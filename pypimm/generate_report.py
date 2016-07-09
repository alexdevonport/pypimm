__author__ = 'alex'
import os, sys, re, subprocess, pandas

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def generate_report(analysis, templatedir):
    """

    :param name:
    :param data:
    :param fits:
    :param properties:
    :return:
    """

    # Get relevant data from the analysis object
    name = analysis.get_name()
    data = analysis.get_raw_data()
    fits = analysis.get_fits()
    properties = analysis.get_results()

    # This is just to get a sense of what pandas output looks like
    fitsFrame = pandas.DataFrame.from_dict(fits).transpose()
    dataFrame = pandas.DataFrame.from_dict(data)
    propertiesFrame = pandas.DataFrame.from_dict(properties)

    writer = pandas.ExcelWriter('RESULTS_'+name+'.xlsx')

    dataFrame.to_excel(writer,sheet_name='fitted data')
    fitsFrame.to_excel(writer,sheet_name='Fit results')
    propertiesFrame.to_excel(writer,sheet_name='calculated properties')

    writer.save()

    # TODO: figure out a way to find out which directory the original
    # TODO: PyPIMM.py file was called from. This should not be hard-coded.

    # get the template file
    template_file = open(os.path.join(templatedir, 'report_template.tex'))
    report_file = open(os.path.join('.',name+'_report.tex'), 'w')

    latexKeyword = re.compile(r'\\KEYWORD\((\w*)\)')

    # This dicitonary will contain the substitutions we'll be making
    # into the template report file.
    substitutions = {}

    # TODO: get name of experiment
    substitutions['experimentName'] = name

    # TODO: get filepaths for summary images

    substitutions['fVersusH'] = r'./'+name+'-f-vs-h.png'
    substitutions['dVersusH'] = r'./'+name+'-d-vs-h.png'

    # TODO: generate lists of images to plot
    fitList = listdir_fullpath(os.path.join('.','sigfits')) # get a list of all filepath strings
    fitList = str(fitList)  # convert the list into a string
    fitList = fitList[1:-1] # remove the '[' and ']' characters at beginning and end of string
    fitList = fitList.replace('\'', '')  # remove single quotes from around items

    substitutions['fitList'] = fitList

    spectList = listdir_fullpath(os.path.join('.','spectra'))
    spectList = str(spectList)
    spectList = spectList[1:-1]
    spectList = spectList.replace('\'', '')
    substitutions['spectList'] = spectList

    template = template_file.read()
    template_subbed = latexKeyword.sub(lambda mo: getValue(mo, substitutions), template)
    report_file.write(template_subbed)
    report_file.close()
    template_file.close()

    # call pdflatex to generate the pdf
    #subprocess.call(['pdflatex', name+'_report.tex'])

    return None

def getValue(mo, subdict):
    kw = mo.group(1)
    try:                 # see if we defined a substitution matching this keyword
        r = subdict[kw]  # if we did, return it
    except KeyError:     # otherwise,
        r = kw           # just leave the keyword
    return r