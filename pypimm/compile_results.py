import pandas as pd
__author__ = 'alex'


def compile_results(rs):
    df = pd.DataFrame()
    usedNames = []
    k = 1
    writer = pd.ExcelWriter('overall-report.xlsx')
    for r in rs:
        name = r['name']
        newname = name
        while newname in usedNames:
            newname = name + ' ({:d})'.format(k)
            k += 1
        usedNames.append(newname)
        df[newname] = pd.Series(r['results'])
    df = df.transpose()
    df.to_excel(writer)
    writer.save()
    return None

#test
def main():
    q = [
            {
                'name':'test set a',
                'results':{'a':1, 'b':2},
            },
            {
                'name':'test set a',
                'results':{'a':3, 'b':4},
            }
        ]
    compile_results(q)

if __name__ == '__main__':
    main()