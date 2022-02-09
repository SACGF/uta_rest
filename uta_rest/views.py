import re

from django.conf import settings
from django.db import connection
from django.http import Http404, JsonResponse
from django.shortcuts import render
from django.views.decorators.cache import cache_page


def index(request):
    return render(request, 'index.html')


@cache_page(settings.DEFAULT_TIMEOUT)
def transcript(request, transcript_version):
    return uta_transcript(request, settings.UTA_DEFAULT_VERSION, transcript_version)


@cache_page(settings.DEFAULT_TIMEOUT)
def uta_transcript(request, uta_version, transcript_version):
    # uta_version needs to look like 'uta_20210129'
    # But connection is a read only public database so don't need to worry too much about SQL injection
    pattern = r"^uta_\d{8}$"
    if not re.match(pattern, uta_version):
        raise ValueError(f"uta_version must be conform to regex: '{pattern}'")

    data = _get_data_from_uta(uta_version, transcript_version)
    if data is None:
        raise Http404(f"transcript '{transcript_version}' not found in UTA '{uta_version}'")

    return JsonResponse(data)


def _get_data_from_uta(uta_version, transcript_version):
    data = {}
    with connection.cursor() as cursor:
        transcript_sql = f"""
        SELECT ac, hgnc, cds_start_i, cds_end_i
        FROM "{uta_version}".transcript transcript 
        WHERE transcript.ac = %s;
        """
        cursor.execute(transcript_sql, [transcript_version])
        if row := cursor.fetchone():
            columns = ["id", "gene_name", "start_codon", "stop_codon"]
            for k, v in zip(columns, row):
                data[k] = v
        else:
            return {}

        data["genome_builds"] = {}

        build_sql = f"""
        SELECT
        string_agg(distinct es.alt_ac::varchar, ',') as contig,
        string_agg(distinct es.alt_strand::varchar, ',') as strand,
        string_agg(exon.start_i::varchar, ',' order by exon.ord) as exon_starts,
        string_agg(exon.end_i::varchar, ',' order by exon.ord) as exon_ends,
        string_agg(exon_aln.cigar, ',' order by exon.ord) as cigars
        from "{uta_version}".transcript transcript
        inner join "{uta_version}".exon_set es on (transcript.ac = es.tx_ac AND alt_aln_method = 'splign')
        inner join "{uta_version}".origin origin on (transcript.origin_id = origin.origin_id)
        Inner join "{uta_version}".exon as exon on (es.exon_set_id = exon.exon_set_id)
        inner join "{uta_version}".exon_aln exon_aln on (exon_aln.alt_exon_id = exon.exon_id)
        WHERE transcript.ac = %s AND es.alt_ac = ANY(%s)
        group by transcript.ac
        """

        for genome_build, contigs in settings.GENOME_BUILD_CONTIGS.items():
            cursor.execute(build_sql, [transcript_version, contigs])
            if row := cursor.fetchone():
                (contig, strand, exon_starts, exon_ends, cigars) = row
                data["genome_builds"][genome_build] = {
                    "contig": contig,
                    "strand": "+" if strand == 1 else "-",
                    "exons": _convert_uta_exons(exon_starts, exon_ends, cigars),
                }
    return data


# Copy/pasted from cdot: https://github.com/SACGF/cdot/blob/main/generate_transcript_data/cdot_json.py
def _convert_uta_exons(exon_starts, exon_ends, cigars):
    # UTA is output sorted in exon order (stranded)
    exon_starts = exon_starts.split(",")
    exon_ends = exon_ends.split(",")
    cigars = cigars.split(",")
    exons = []
    ex_ord = 0
    ex_transcript_start = 1  # transcript coords are 1 based
    for ex_start, ex_end, ex_cigar in zip(exon_starts, exon_ends, cigars):
        gap, exon_length = _cigar_to_gap_and_length(ex_cigar)
        ex_transcript_end = ex_transcript_start + exon_length - 1
        exons.append((int(ex_start), int(ex_end), ex_ord, ex_transcript_start, ex_transcript_end, gap))
        ex_transcript_start += exon_length
        ex_ord += 1

    exons.sort(key=lambda e: e[0])  # Genomic order
    return exons

def _cigar_to_gap_and_length(cigar):
    """
            gap = 'M196 I1 M61 I1 M181'
            CIGAR = '194=1D60=1D184='
    """

    # This has to/from sequences inverted, so insertion is a deletion
    OP_CONVERSION = {
        "=": "M",
        "D": "I",
        "I": "D",
        "X": "=",  # TODO: This is probably wrong! check if GTF gap has mismatch?
    }

    cigar_pattern = re.compile(r"(\d+)([" + "".join(OP_CONVERSION.keys()) + "])")
    gap_ops = []
    exon_length = 0
    for (length_str, cigar_code) in cigar_pattern.findall(cigar):
        exon_length += int(length_str)  # This shouldn't include one of the indels?
        gap_ops.append(OP_CONVERSION[cigar_code] + length_str)

    gap = " ".join(gap_ops)
    if len(gap_ops) == 1:
        gap_op = gap_ops[0]
        if gap_op.startswith("M"):
            gap = None

    return gap, exon_length
