import argparse, glob, os
from Bio import SeqIO

__desc__ = """TATA search looks for the TATA sequence in a FASTA file
              if the TATA sequence is within a designated number of
              bases from a start codon and if the ORF that follows it
              is of a designated length. TATA search produces an html
              file that contains a list of all found sequences, their
              lengths and the protein sequence of the ORF as a report."""


def DNA_rev_comp(sekvenca):
    for letter,repl in zip('ATGC10','1A0GTC'):
        sekvenca=sekvenca.replace(letter,repl)
    return sekvenca[::-1]


def parse_arguments():
    """
    Parse command-line arguments using argparse module.
    :return: Filepath to the input FASTA file
    """
    parser = argparse.ArgumentParser(description=__desc__)

    parser.add_argument('fasta',
                        type=str,
                        default='',
                        help='FASTA file (.fasta is optional)')
    parser.add_argument('-t', '--tata',
                        type=int,
                        default=25,
                        help='TATA Distance from START (default 25). By default TATA search looks for distances of +/- 10% of the given value')
    parser.add_argument('-o', '--orf',
                        type=int,
                        default=0,
                        help='Minimum ORF length in bp (default 0)')
    parser.add_argument('-a','--allfasta',
                        action='store_true',
                        default=False,
                        help='Find and read all FASTA files (default False)')
    args = vars(parser.parse_args())

    return args['fasta'], args['tata'], args['orf'],args['allfasta']


def find_all(iterable, strings):
    """
    Searches for all indices of "string" in the "iterable".
    :param iterable: A listable variable.
    :param strings: A string is parsed as a single input, while a list as multiple inputs
    :return: Indices of string in the iterable
    """
    lokacije = []
    start = 0
    if type(strings) == str:
        strings = [strings]

    for string in strings:
        while True:
            lokacija = iterable[start:].find(string)
            if lokacija >= 0:
                lokacije.append(lokacija + start)
                start += lokacija + 1
            else:
                break
    return lokacije


def main():
    # inputs
    files, tata, orf, allfasta = parse_arguments()

    if not files or allfasta:
        files=[]
        # os.chdir(os.getcwd())
        for file in glob.glob('*.fasta'):
            files.append(file)
    else:
        files=files.split(',')

    duzina_okvira = orf
    lokacija_TATA = tata

    # input prep
    for i,f in enumerate(files):
        if '.fasta' not in f:
            files[i] = f+'.fasta'

    # output prep
    outputi = []

    # script
    for file in files:
        html_name=file+'.html'
        i = SeqIO.read(file, "fasta")
        for loc1 in find_all(i.seq, 'ATG'):
            for loc2 in find_all(i.seq, 'TATA'):
                if loc1 - lokacija_TATA * 1.1 < loc2 < loc1 - lokacija_TATA * 0.9:
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
                        outputi.append('</td><td>'.join(str(deo) for deo in
                                                 [loc2, loc2 + len(sekvenca), len(sekvenca), i.seq[loc2:loc2 + 8]+
                                                  '... ' + str(loc1 - loc2 - 8) + ' ...'+ sekvenca,sekvenca.translate()]))

        # output
        message = """
        <html>
        <head>
        <h3>Outputi za {file}</h3>
        <style>
            table {{
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }}
        td, th {{
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }}
        tr:nth-child(even) {{
            background-color: #dddddd;
        }}
        </style>
        </head>
        <body>
        <table>
          <tr>
            <th>From</th>
            <th>To</th>
            <th>Length</th>
            <th>DNA sequence</th>
            <th>AA sequence</th>
          </tr>""".format(file=file) + ''.join(
            ["<tr><td>" + out + "</td></tr>" for out in outputi]) + "</table></body></html>"

        with open(html_name, 'w') as prikaz:
            prikaz.write(message)


if __name__ == '__main__':
    main()
