import requests
from bs4 import BeautifulSoup
import sys
num_genomes = 30
if len(sys.argv)>1:
    num_genomes = int(sys.argv[1])
top = []
with open("datarank.txt",'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.split(',')
        name = split[0]
        top.append(name)
top = top[0:num_genomes]
from collections import OrderedDict
url = 'https://hgdownload.soe.ucsc.edu/hubs/birds/index.html'
r = requests.get(url)
html_content = r.text
soup = BeautifulSoup(html_content, 'lxml')
links = []
for link in soup.findAll('a'):
    if "hgdownload" in str(link):
        links.append(link.get('href'))
links = links[7::]
theLinks = []
for l in links:
    s = l.split('/')
    name = s[8]
    if name in top:
        print(name)
        theLinks.append(l)
f = open("fasta_links_small.txt",'w')
for l in theLinks:
    #print(l)
    req = requests.get(l)
    html_content2 = req.text
    soup = BeautifulSoup(html_content2, 'lxml')
    for link in soup.findAll('a'):
        if "fa.gz" in str(link):
            ls = link.get('href')
            print(str(l)+str(ls))
            f.write(str(l)+str(ls)+'\n')


