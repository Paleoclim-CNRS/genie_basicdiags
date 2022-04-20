#!/usr/bin/env python

def dopdf(savedfilesdata, prefixdata, filecountdata, namedata, sentence_to_be_replaced_data):

    import os
    from shutil import copyfile
    import re

    # copy template
    copyfile('source/empty_v2.tex', 'summary.tex')

    # ===================== TITLE =====================

    documentname = namedata + '_summary.pdf'
    titlename = namedata
 
    # read in the file
    with open('summary.tex', 'r') as file :
        filedata = file.read()
    # replace string
    filedata = filedata.replace('titlestring', titlename)
    filedata = filedata.replace('sentence_to_be_replaced', sentence_to_be_replaced_data)

    # write the file out again
    with open('summary.tex', 'w') as file:
        file.write(filedata)

    # ===================== FIGURES =====================

    counter = 0
    for filecnt in np.arange(filecountdata)+1:
        pdffile = prefixdata + str(filecnt)
        # copy template
        with open('source/figurestatement.tex', 'r') as figurefile :
            figurefiledata = figurefile.read()
        # replace string
        figurefiledata = figurefiledata.replace('figstring', pdffile)
        figurefiledata = figurefiledata.replace('captionstring', re.sub('_', '\_',savedfilesdata[counter][:-4]))
        figurefiledata = figurefiledata.replace('labelstring', re.sub('_', '\_', savedfilesdata[counter][:-4]))
        # add created text into document
        with open('summary.tex', 'a') as texfile:
            texfile.writelines(figurefiledata)
        counter += 1

    # ===================== END OF FILE =====================

    endstatement = '\end{document}'
    with open('summary.tex', 'a') as texfile:
        texfile.writelines(endstatement)

    # ===================== END OF FILE =====================

    os.system("pdflatex summary.tex")
    os.system("pdflatex summary.tex")
    os.remove("summary.out")
    os.remove("summary.aux")
    os.remove("summary.lof")
    os.remove("summary.log")
    os.system("rm -rf file*")
    os.system("rm -rf texput*")
    os.system("mv summary.pdf " + documentname)
    os.remove("summary.tex")
