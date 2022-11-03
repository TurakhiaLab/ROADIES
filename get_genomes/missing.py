import requests
from bs4 import BeautifulSoup
import sys

need = []
with open('need.txt','r') as f:
    lines=f.readlines()
    for l in lines:
        if len(l)> 0:
            need.append(l.strip())
print(need)
fastas = {}
with open('species_map.txt','r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        sn = split[1]
        num = split[0]
        print(sn)
        if sn in need:
            if sn not in fastas:
                fastas[num] = sn
print(fastas)
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
    if name in fastas:
        print(name)
        theLinks.append(l)
f = open("links48.txt",'w')
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



