import argparse
from Bio import SeqIO

__desc__ = """TATA search looks for the TATA sequence in a FASTA file
              if the TATA sequence is within a designated number of
              bases from a start codon and if the ORF that follows it
              is of a designated length. TATA search produces an html
              file that contains a list of all found sequences, their
              lengths and the protein sequence of the ORF as a report."""


def parse_arguments():
    """
    Parse command-line arguments using argparse module.
    :return: Filepath to the input FASTA file
    """
    parser = argparse.ArgumentParser(description=__desc__)
    parser.add_argument('-f', '--fasta',
                        type=str,
                        help='FASTA file (.fasta is optional)',
                        required=True)
    parser.add_argument('-t', '--tata',
                        type=int,
                        help='TATA Distance from START')
    parser.add_argument('-o', '--orf',
                        type=int,
                        help='Minimum ORF length in bp')

    args = vars(parser.parse_args())

    return args['fasta'], args['tata'], args['orf']


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


def main(duzina_okvira=0, lokacija_TATA=30):
    # inputs
    file, tata, orf = parse_arguments()

    if orf:
        duzina_okvira = orf
    if tata:
        lokacija_TATA = tata

    # input prep
    if '.fasta' not in file:
        file += '.fasta'

    # output prep
    outputi = []

    # script
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
                    outputi.append('\t'.join(str(deo) for deo in
                                             [loc2, loc2 + len(sekvenca), len(sekvenca), i.seq[loc2:loc2 + 8],
                                              '... ' + str(loc1 - loc2 - 8) + ' ...', sekvenca, sekvenca.translate()]))

    message = "<html><head><h3>Outputi</h3></head><body>" + ''.join(
        ["<p>" + out + "</p>" for out in outputi]) + "</body></html>"

    with open('prikaz.html', 'w') as prikaz:
        prikaz.write(message)


if __name__ == '__main__':
    main()
