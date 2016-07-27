import pandas as pd
__author__ = 'alex'


def compile_results(rs):
    writer = pd.ExcelWriter('overall-report.xlsx')
    for r in rs:
        name = r['name']
        try:
            name = name[:29]
        except IndexError:
            pass
        result = r['results']
        res_frame = pd.DataFrame.from_dict(result)
        res_frame.to_excel(writer, sheet_name=name)
        writer.save()