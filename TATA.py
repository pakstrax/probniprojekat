from Bio import SeqIO
import webbrowser
from os import getcwd


def find_all(iterable, strings):
    lokacije = []
    start = 0
    if type(strings)==str:
        strings=[strings]

    for string in strings:
        while True:
            lokacija = iterable[start:].find(string)
            if lokacija >= 0:
                lokacije.append(lokacija + start)
                start += lokacija + 1
            else:
                break
    return lokacije


# inputi
duzina_okvira=30
lokacija_TATA=25
i = SeqIO.read("example_human_reference.fasta", "fasta")


outputi=[]

for loc1 in find_all(i.seq, 'ATG'):
    for loc2 in find_all(i.seq, 'TATA'):
        if loc1-lokacija_TATA*1.1< loc2 < loc1-lokacija_TATA*0.9:
            prekidac = True
            start = loc1
            sekvenca = ""
            triplet = ""
            while triplet not in ['TAG', 'TAA', 'TGA']:
                triplet = i.seq[start:start + 3]
                sekvenca += triplet
                start += 3
                if start > len(i.seq):
                    prekidac = False
                    break
            if prekidac and (len(sekvenca) > duzina_okvira):
                outputi.append('\t'.join([loc2,loc2+len(sekvenca), len(sekvenca), i.seq[loc2:loc2+8], '... '+loc1 - loc2-8+' ...',sekvenca, sekvenca.translate()]))

message="<html><head><h3>Outputi</h3></head><body>"+''.join(["<p>"+out+"</p>" for out in outputi])+"</body></html>"

with open('prikaz.html','w') as prikaz:
    prikaz.write(message)

filename=getcwd()+'prikaz.html'
webbrowser.open_new_tab(filename)
