import os
import re
import tempfile
import json
from operator import __or__

from django.db import connection
from django.db.models import Q
from django.http import JsonResponse, HttpResponse
from django.views.decorators.gzip import gzip_page

from .models import Variant
#########################################################
from django.views.decorators.http import require_http_methods
from ga4gh import variant_service_pb2 as v_s
from ga4gh import variants_pb2 as vrs
from ga4gh import metadata_service_pb2 as m_s
from ga4gh import metadata_pb2 as meta
import google.protobuf.json_format as json_format

#########################################################

@gzip_page
def index(request):
    order_by = request.GET.get('order_by')
    direction = request.GET.get('direction')
    page_size = int(request.GET.get('page_size', '0'))
    page_num = int(request.GET.get('page_num', '0'))
    search_term = request.GET.get('search_term')
    format = request.GET.get('format')
    include = request.GET.getlist('include')
    exclude = request.GET.getlist('exclude')
    filters = request.GET.getlist('filter')
    filter_values = request.GET.getlist('filterValue')
    column = request.GET.getlist('column')

    query = Variant.objects

    if format == 'csv':
        quotes = '\''
    else:
        quotes = ''

    if include or exclude:
        query = apply_sources(query, include, exclude)

    if filters:
        query = apply_filters(query, filter_values, filters, quotes=quotes)

    if search_term:
        query = apply_search(query, search_term, quotes=quotes)

    if order_by:
        query = apply_order(query, order_by, direction)

    if format == 'csv':

        cursor = connection.cursor()
        with tempfile.NamedTemporaryFile() as f:
            os.chmod(f.name, 0606)
            cursor.execute("COPY ({}) TO '{}' WITH DELIMITER ',' CSV HEADER".format(query.query, f.name))

            response = HttpResponse(f.read(), content_type='text/csv')
            response['Content-Disposition'] = 'attachment;filename="variants.csv"'
            return response

    elif format == 'json':

        count = query.count()

        if search_term:
            # Number of synonym matches = total matches minus matches on "normal" columns
            synonyms = count - apply_search(query, search_term, search_column='fts_standard').count()
        else:
            synonyms = 0

        query = select_page(query, page_size, page_num)

        # call list() now to evaluate the query
        response = JsonResponse({'count': count, 'synonyms': synonyms, 'data': list(query.values(*column))})
        response['Access-Control-Allow-Origin'] = '*'
        return response


def apply_sources(query, include, exclude):
    # if there are multiple sources given then OR them:
    # the row must match in at least one column
    include_list = (Q(**{column: True}) for column in include)
    exclude_dict = {exclusion: False for exclusion in exclude}

    return query.filter(reduce(__or__, include_list)).filter(**exclude_dict)


def apply_filters(query, filterValues, filters, quotes=''):
    # if there are multiple filters the row must match all the filters
    for column, value in zip(filters, filterValues):
        if column == 'id':
            query = query.filter(**{column: value})
        else:
            query = query.extra(
                where=["\"{0}\" LIKE %s".format(column)],
                params=["{0}{1}%{0}".format(quotes, value)]
            )
    return query


def apply_search(query, search_term, search_column='fts_document', quotes=''):
    # search using the tsvector column which represents our document made of all the columns
    where_clause = "variant.{} @@ to_tsquery('simple', %s)".format(search_column)
    parameter = quotes + sanitise_term(search_term) + quotes
    return query.extra(
        where=[where_clause],
        params=[parameter]
    )


def apply_order(query, order_by, direction):
    # special case for HGVS columns
    if order_by in ('HGVS_cDNA', 'HGVS_Protein'):
        order_by = 'Genomic_Coordinate_hg38'
    if direction == 'descending':
        order_by = '-' + order_by
    return query.order_by(order_by, 'Pathogenicity_default')


def select_page(query, page_size, page_num):
    if page_size:
        start = page_size * page_num
        end = start + page_size
        return query[start:end]
    return query


def autocomplete(request):
    term = request.GET.get('term')
    limit = int(request.GET.get('limit', 10))

    cursor = connection.cursor()

    cursor.execute(
        """SELECT word FROM words
        WHERE word LIKE %s
        AND char_length(word) >= 3
        ORDER BY word""",
        ["%s%%" % term])

    rows = cursor.fetchall()

    response = JsonResponse({'suggestions': rows[:limit]})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def sanitise_term(term):
    # Escape all non alphanumeric characters
    term = re.escape(term)
    # Enable prefix search
    term += ":*"
    return term

########################### START WORK ####################################
###########################################################################
###########################################################################

############## .../variants/search method ############################
@require_http_methods(["POST"])
def variant_search(request):
    conditional = validate_request(request)
    if conditional :
        return conditional
    else:
        req_dict = json.loads(request.body)
        variant_set_id = req_dict.get('variantSetId')
        reference_name = req_dict.get('referenceName')
        start = req_dict.get('start')
        end = req_dict.get('end')
        page_size = req_dict.get('pageSize', 3)
        page_token = req_dict.get('pageToken', "0")

    x, referenceCord = variant_set_id.split("-")
    if referenceCord not in SetIds :
        return JsonResponse({"Error messeage" : "Not a supported variantSetId"})

    response0 = v_s.SearchVariantsResponse()
    filt = str(reference_name)+":"
    DbResp = Variant.objects
    DbResp = select_gen_coor(referenceCord, DbResp, filt)

    ret_data = ga4gh_brca_page(DbResp, int(page_size), int(page_token))

    ga_vars = [brca_to_ga4gh(i,referenceCord) for i in ret_data.values()]
    if len(ga_vars) > page_size:
        ga_vars.pop()
        page_token = str(1 + int(page_token))
    else:
        response0.next_page_token = ""
    response0.variants.extend(ga_vars)
    response0.next_page_token = page_token
    resp = json_format._MessageToJsonObject(response0, True)
    return JsonResponse(resp)
############### .../variants/search methond END ##############################

############### Auxiliary function for quering database on a given Gen-Coordinate ###############
def select_gen_coor(referenceCordinate, DbResp, filt, Start = None, End = None):
    if referenceCordinate == SetIds[0]:
        FilteredResponse = DbResp.filter(Genomic_Coordinate_hg36__startswith=filt)
        if Start != None and End != None:
            FilteredResponse = FilteredResponse.filter(Hg36_End__lt=End, Hg36_Start__gt=val1)
        else:
            return FilteredResponse
    elif referenceCordinate == SetIds[1]:
        FilteredResponse = DbResp.filter(Genomic_Coordinate_hg37__startswith=filt)
        if Start != None and End != None:
            FilteredResponse = FilteredResponse.filter(Hg37_End__lt=End, Hg37_Start__gt=val1)
        else:
            return FilteredResponse
    elif referenceCordinate == SetIds[2]:
        FilteredResponse = DbResp.filter(Genomic_Coordinate_hg38__startswith=filt)
        if Start != None and End != None:
            FilteredResponse = FilteredResponse.filter(Hg38_End__lt=End, Hg38_Start__gt=val1)
        else:
            return FilteredResponse
############################## Gen-Coordinate function END ##############################

############### Paging Auxiliary function ###############
def ga4gh_brca_page(query, page_size, page_token):
    start = page_size * page_token
    end = start + page_size+1
    return query[start:end]
################## Paging function END ##################

############### Auxiliary function to obtain offset for database query ########
def get_offset(start, end, VarLEn = None):
    if VarLEn:
        start = start - 1
        end = start + len(VarLEn)
    else:
        start = start + 1
        end = end + 1
    return start, end
###################### Offset function END ##############################


### Function that translates elements in BRCA-database to GA4GH format #######
def brca_to_ga4gh(brca_variant, genRefer):
    var_resp = vrs.Variant()
    for j in brca_variant:
        if brca_variant[j] == "-" or (brca_variant[j] == ""):
            continue
        if j == "Genomic_Coordinate_"+genRefer:
            var_resp.reference_name, start, bases = brca_variant[j].split(':')
            var_resp.reference_bases, alternbases = bases.split(">")
            for i in range(len(alternbases)):
                var_resp.alternate_bases.append(alternbases[i])
            var_resp.start = int(start)
            var_resp.end = var_resp.start+len(alternbases)
            continue
        if j == "id":
            var_resp.id = genRefer+"-"+str(brca_variant['id'])
            var_resp.variant_set_id = name+"-"+genRefer

            var_resp.created = 0
            var_resp.updated = 0
            continue
        if j == "Synonyms":
            Names = [i for i in str(brca_variant[j]).split(",")]
            for i in Names:
                var_resp.names.append(i)
        else:
            var_resp.info[str(j)].append(brca_variant[j])
    return var_resp
################## Translator function END ##############################

#### Auxiliary function which validates requests ######
def validate_request(request):
    if not request.body:
        return JsonResponse(ErrorMessages['emptyBody'])
    else:
        request_dict = json.loads(request.body)
        if not request_dict.get("variantSetId"):
            return JsonResponse(ErrorMessages['VariantSetId'])
        elif not request_dict.get('referenceName'):
            return JsonResponse(ErrorMessages['referenceName'])
        elif not request_dict.get('start')  :
            return JsonResponse(ErrorMessages['start'])
        elif not request_dict.get('end') :
            return JsonResponse(ErrorMessages['end'])
        else:
            return None

def validate_varsetreq(request):
    if not request.body:
        return JsonResponse(ErrorMessages['emptyBody'])
    else:
        request_dict = json.loads(request.body)
        if not request_dict.get("datasetId"):
            return JsonResponse(ErrorMessages['datasetId'])
        else:
            return None
########## Auxiliary function END ############

######### .../variants/<variant id> method ########
@require_http_methods(["GET"])
def get_var_by_id(request, variant_id):
    if not variant_id:
        return JsonResponse(ErrorMessages['variantId'])
    else:
        gen_coor_and_id = variant_id
        gen_coor, v_id = gen_coor_and_id.split("-")
        if gen_coor in SetIds:
            DbResp = Variant.objects.values()
            resp1 = DbResp.get(id=int(v_id))
            Var_resp = brca_to_ga4gh(resp1, gen_coor)
            resp = json_format._MessageToJsonObject(Var_resp, True)
            return JsonResponse(resp)
        else:
            return JsonResponse(ErrorMessages["datasetId"])
######### .../variants/<variant id> method END ########

######### .../variantsets/search method ########
@require_http_methods(["POST"])
def get_variantSet(request):
    condit = validate_varsetreq(request)
    if condit:
        return condit
    else:
        req_dict = json.loads(request.body)
        dataset_id = req_dict.get('datasetId')
        page_size = req_dict.get('pageSize', 3)
        page_token = req_dict.get('pageToken', '0')
        if dataset_id != "brca-exchange":
            return JsonResponse(ErrorMessages["datasetId"])
    if page_token is None:
         page_token = "0"

    response1 = v_s.SearchVariantSetsResponse()
    response1.next_page_token = page_token
    for i in range(len(SetIds)):
        response = vrs.VariantSet()
        response.id = datasetId+"-"+SetIds[i]
        response.name = SetName+"-"+SetIds[i]
        response.dataset_id =  name
        response.reference_set_id = referenceSetId+"-"+SetIds[i]
        brca_meta(response.metadata, dataset_id)
        response1.variant_sets.extend([response])
    resp = json_format._MessageToJsonObject(response1, True)
    return JsonResponse(resp)
######### .../variantsets/search method END ########

### Auxiliary function, generates metadata fields ######
def brca_meta(Metadata, dataset_id):
    var_resp = vrs.VariantSetMetadata()
    for j in Variant._meta.get_all_field_names():
        var_resp.key = str(j)
        var_resp.value = "-"
        var_resp.id = datasetId+"-"+str(j)
        var_resp.type = Variant._meta.get_field(str(j)).get_internal_type()
        var_resp.number = "-"
        var_resp.description = "refer to ->"+str(j)+ " in https://github.com/BD2KGenomics/brca-website/blob/master/content/help_research.md"
        Metadata.extend([var_resp])
    return Metadata
### Metadata generator END ###############

#### .../variantsets/<set id> method ########
@require_http_methods(["GET"])
def get_varset_by_id(request, variantSetId):
    if not variantSetId :
        return JsonResponse(ErrorMessages["variantSetId"])
    dataset, Id = variantSetId.split("-")
    if Id in SetIds and dataset == "brca":
        response = vrs.VariantSet()
        response.id = dataset + "-" + Id
        response.name = SetName + "-" + Id
        response.dataset_id = name
        response.reference_set_id = referenceSetId + "-" + Id
        brca_meta(response.metadata, Id)
        resp = json_format._MessageToJsonObject(response, True)
        return JsonResponse(resp)
    else:
        return JsonResponse({"Invalid Set Id": variantSetId})
#### .../variantsets/<set id> method END ########

#### .../datasets/search method request handler ######
@require_http_methods(["POST"])
def search_datasets(request):
    req_dict = json.loads(request.body)
    page_size = req_dict.get("pageSize", 1)
    page_token = req_dict.get("nextPageToken", "0")
    sr_dta_set_resp = m_s.SearchDatasetsResponse()
    dta_resp = meta.Dataset()
    dta_resp.name = SetName
    dta_resp.id = name
    dta_resp.info[SetName].append("This set contains variants as stored and mantained by the brca-exchange project")
    dta_resp.description = "Variants observed in brca-exchange project"
    sr_dta_set_resp.datasets.extend([dta_resp])
    sr_dta_set_resp.next_page_token = page_token
    return JsonResponse(json_format._MessageToJsonObject(sr_dta_set_resp, True))
#### .../datasets/search method END #########

 #### Error URL catcher methods ######
@require_http_methods(["GET", "POST"])
def varsetId_empty_catcher(request):
    return JsonResponse(ErrorMessages["emptyBody"])

@require_http_methods(["GET", "POST"])
def empty_varId_catcher(request):
    return JsonResponse(ErrorMessages["emptyBody"])
################# END #####################

############################## GLOABAL VARIABLES ################################
ErrorMessages = {'emptyBody' :{'error code': 400, 'message' : 'invalid request empty request'},
                'VariantSetId' : {'error code': 400, 'message': 'invalid request no variant_set_id'},
                 'referenceName': {'error code': 400, 'message': 'invalid reques no reference_name'},
                 'start': {'error code' : 400, 'message': 'invalid request no start'},
                 'end' : {'error code' :400, 'message': 'invalid request no end'},
                 'datasetId': {'error code' : 400, 'message': 'invalid request no dataset_id'},
                 'variantId': {'error code' : 400, 'message': 'invalid request no variant_id'},
                 'variantSetId': {'error code': 400, 'message': 'invalid request no variant_set_id'}}
SetName = "brca-exchange-variants"
datasetId = "brca"
referenceSetId = "Genomic-Coordinate"
name = "brca-exchange"
SetIds = ["hg36", "hg37", "hg38"]
#################################################################################

# Need to implement function that filters elements by increasing and decreasing integer values
