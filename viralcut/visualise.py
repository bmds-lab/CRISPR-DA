from pathlib import Path

from . import data
from . import utils
from .collection import ViralCutCollection
from ete3.parser.newick import write_newick

def generate_itol_tree(collection: ViralCutCollection):
    outputDir = Path('/mnt/ssd1/carl/ViralCut/outputs') # TODO: Replace this with user defined dir
    outputDir.mkdir(exist_ok=True)

    with open(outputDir / 'tree.nwk', 'w') as outFile:
        outFile.write(write_newick(collection._ncbi_tree))

    taxIdToName = {}
    for accession, name in zip(collection.accessions, data.get_name_from_accessions(collection.accessions)):
        taxId = str(collection.accession_to_tax_id[accession])
        if taxId in taxIdToName:
            assert(name == taxIdToName[taxId])
        else:
            taxIdToName[taxId] =  name
    leaves = collection._ncbi_tree.get_leaves()

    for guide in [y for y, x in collection.guides.items() if x.assembly_scores != {}]:
        guideDir = outputDir / guide
        guideDir.mkdir(exist_ok=True)

        with open(guideDir / 'labels.txt', 'w') as outFile:
            outFile.write('LABELS\nSEPARATOR COMMA\nDATA\n')
            for leaf in leaves:
                outFile.write(f'{leaf.taxid},{taxIdToName[str(leaf.taxid)]}\n')

        with open(guideDir / 'popups.txt', 'w') as outFile:
            outFile.write('POPUP_INFO\nSEPARATOR COMMA\nDATA\n')
            for leaf in leaves:
                try:
                    name = taxIdToName[str(leaf.taxid)]
                except:
                    name = f'Group of {len(leaf.children)} nodes'
                outFile.write(f'{leaf.taxid},' +  f'{name},' +
                            f'<p>MIT Score: {leaf.score[guide]['mit']}</p>' +
                            f'<p>CFD Score: {leaf.score[guide]['cfd']}</p>' +
                            f'<p>Total Off-target sties: {leaf.score[guide]['total_sites']}</p>' +
                            f'<p>Unique Off-target sties: {leaf.score[guide]['unique_sites']}</p>\n')
            for node in collection._ncbi_tree.traverse("preorder"):
                if len(node.children) > 0:
                    firstChild = node.children[0]
                    while len(firstChild.children) > 0:
                        firstChild = firstChild.children[0]
                    lastChild = node.children[-1]
                    while len(lastChild.children) > 0:
                        lastChild = lastChild.children[-1]
                    try:
                        name = taxIdToName[str(node.taxid)]
                    except:
                        name = f'Group of {len(node.children)} nodes'
                    outFile.write(f'{firstChild.taxid}|{lastChild.taxid},' +  f'{name},' +
                            f'<p>MIT Score: {node.score[guide]['mit']}</p>' +
                            f'<p>CFD Score: {node.score[guide]['cfd']}</p>' +
                            f'<p>Total Off-target sties: {node.score[guide]['total_sites']}</p>' +
                            f'<p>Unique Off-target sties{node.score[guide]['unique_sites']}</p>\n')

        with open(guideDir / 'MIT.txt', 'w') as outFile:
            outFile.write('DATASET_COLORSTRIP\n' +
                        'SEPARATOR COMMA\n' + 
                        'DATASET_LABEL,MIT\n' +
                        'COLOR,#ff0000\n' + 
                        'LEGEND_TITLE,Off-target scores\n'+
                        'LEGEND_SHAPES,1,1,1,1\n' +
                        'LEGEND_COLORS,#e06666,#f6b26b,#ffe599,#93c47d\n' + 
                        'LEGEND_LABELS,0-25,25-50,50-75,75-100\n' +
                        'DATA\n'
                        )
            for leaf in leaves:
                outFile.write(f'{leaf.taxid},{utils.colour_code(leaf.score[guide]['mit'])},{leaf.score[guide]['mit']}\n')
            for node in collection._ncbi_tree.traverse("preorder"):
                if len(node.children) > 0:
                    firstChild = node.children[0]
                    while len(firstChild.children) > 0:
                        firstChild = firstChild.children[0]
                    lastChild = node.children[-1]
                    while len(lastChild.children) > 0:
                        lastChild = lastChild.children[-1]
                    outFile.write(f'{firstChild.taxid}|{lastChild.taxid},{utils.colour_code(node.score[guide]['mit'])},{node.score[guide]['mit']}\n')

        with open(guideDir / 'CFD.txt', 'w') as outFile:
            outFile.write('DATASET_COLORSTRIP\n' +
                    'SEPARATOR COMMA\n' + 
                    'DATASET_LABEL,CFD\n' +
                    'COLOR,#ff0000\n' + 
                    'LEGEND_TITLE,Off-target scores\n'+
                    'LEGEND_SHAPES,1,1,1,1\n' +
                    'LEGEND_COLORS,#e06666,#f6b26b,#ffe599,#93c47d\n' + 
                    'LEGEND_LABELS,0-25,25-50,50-75,75-100\n' +
                    'DATA\n' 
                    )
            for leaf in leaves:
                outFile.write(f'{leaf.taxid},{utils.colour_code(leaf.score[guide]['cfd'])},{leaf.score[guide]['cfd']}\n')
            for node in collection._ncbi_tree.traverse("preorder"):
                if len(node.children) > 0:
                    firstChild = node.children[0]
                    while len(firstChild.children) > 0:
                        firstChild = firstChild.children[0]
                    lastChild = node.children[-1]
                    while len(lastChild.children) > 0:
                        lastChild = lastChild.children[-1]
                    outFile.write(f'{firstChild.taxid}|{lastChild.taxid},{utils.colour_code(node.score[guide]['cfd'])},{node.score[guide]['cfd']}\n')
