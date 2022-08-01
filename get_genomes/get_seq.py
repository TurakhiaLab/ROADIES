import requests
from bs4 import BeautifulSoup
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
f = open("fasta_links.txt",'w')
for l in links:
    #print(l)
    req = requests.get(l)
    html_content2 = req.text
    soup = BeautifulSoup(html_content2, 'lxml')
    for link in soup.findAll('a'):
        if "fa.gz" in str(link):
            ls = link.get('href')
            print(str(l)+str(ls))
            f.write(str(l)+str(ls)+'\n')


