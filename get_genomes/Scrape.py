import requests
from bs4 import BeautifulSoup
from collections import OrderedDict
url = 'https://hgdownload.soe.ucsc.edu/hubs/birds/index.html'
r = requests.get(url)
html_content = r.text
soup = BeautifulSoup(html_content, 'lxml')
links = []
for link in soup.findAll('a'):
    if "assembly" in str(link).split('/'):
        links.append(link.get('href'))
print(links)
names = []
for link in links:
    split = (link.split('/'))
    name = split[4]
    url = link
    r = requests.get(url)
    html_content = r.text
    soup = BeautifulSoup(html_content)
    td = soup.findAll('td')
    for d in range(len(td)):
        if 'Contig N50' in str(td[d]):
            s = str(td[d+1]).split('>')
            s2 = s[1].split('<')
            num = int(s2[0].replace(',',''))
            names.append((name,num,str(link)))

sorted_by_value = sorted(names,key=lambda tup:tup[1],reverse=True)
with open('datarank.txt','w') as f:
    for d in sorted_by_value:
        f.write(d[0]+',' + str(d[1])+ ',' + d[2]+'\n')

