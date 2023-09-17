import sys
import requests
from requests.adapters import HTTPAdapter, Retry
import re



species = sys.argv[1]

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Cgo_id%2Cgene_oln%2Cgene_primary%2Cgene_synonym%2Cgene_orf%2Cxref_geneid%2Cxref_genetree&format=tsv&query=%28{0}%29'.format(species)


all_fastas = requests.get(url).text
print(all_fastas)

