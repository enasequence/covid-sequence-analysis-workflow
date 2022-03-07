#!/usr/bin/env python3

import sys, os
import xmltodict
import json, requests

def get_sample_xml(run_id):
    adv_search_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22{}%22&fields=sample_accession&format=json".format(run_id)
    response = requests.get(adv_search_url)
    data = json.loads(response.content)
    sample_acc = data[0]['sample_accession']

    portal_xml_url = "https://www.ebi.ac.uk/ena/browser/api/xml/{}".format(sample_acc)
    response = requests.get(portal_xml_url)
    return response.text

def get_attr(attrs, attr_name):
    try:
        return attrs[attr_name]
    except KeyError:
        return 'not provided'

def get_sample_attributes(sample_xml):
    sample_dict = xmltodict.parse(sample_xml)
    sample_acc = sample_dict['SAMPLE_SET']['SAMPLE']['@accession']

    sample_attributes_raw = {}
    if type(sample_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']) is list:
        sample_attributes_raw = {attr['TAG']:attr['VALUE'] for attr in sample_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']}
    else:
        attr = sample_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
        sample_attributes_raw[attr['TAG']] = attr['VALUE']

    # normalise casing and spacing of metadata fields
    sample_attributes = {k.lower().replace(' ', '_'): v for k, v in sample_attributes_raw.items()}

    attributes = {}
    attributes['isolate'] = get_attr(sample_attributes, 'isolate')
    attributes['collection_date'] = get_attr(sample_attributes, 'collection_date')

    if sample_acc[:3] == 'SRS':
        attributes['location'] = get_attr(sample_attributes, 'geo_loc_name')
    else:
        attributes['location'] = get_attr(sample_attributes, 'geographic_location_(country_and/or_sea)')
        try:
            region = sample_attributes['geographic_location_(region_and_locality)']
            if not region[:4] == 'not ':
                attributes['location'] += f":{region}"
        except KeyError:
            pass

    return attributes

def replace_header(run_id):
    sample_xml = get_sample_xml(run_id)
    sample_attributes = get_sample_attributes(sample_xml)

    # format new header
    new_header = ">CDP:{} {}".format( run_id, '|'.join([
        sample_attributes['isolate'], sample_attributes['location'],
        sample_attributes['collection_date'],
        "SARS-CoV-2 complete genome produced by the VEO task force from public INSDC data",
        "covid_sequence_analysis_workflow"
    ]))
    return new_header

with open(sys.argv[1], 'r') as fasta:
    for line in fasta:
        line = line.strip()
        try:
            if line[0] == '>':
                print(replace_header(line[1:]))
            else:
                print(line)
        except IndexError:
            pass
