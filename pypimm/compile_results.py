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

def compile_damping(rs):
    df = pd.DataFrame()
    usedNames = []
    k = 1
    writer = pd.ExcelWriter('damping-report.xlsx')
    for r in rs:
        name = r['name']
        newname = name
        while newname in usedNames:
            newname = name + ' ({:d})'.format(k)
            k += 1
        usedNames.append(newname)
        df[newname + ' field'] = pd.Series(r['results']['hs'])
        df[newname + ' damping'] = pd.Series(r['results']['Damping'])
    df.to_excel(writer)
    writer.save()

    return None

#test
def main():
    q = [
            {
                'name':'test set a',
                'results':{'a':1, 'b':2, 'damping':[1,2,3], 'hs':[1,2,3]},
            },
            {
                'name':'test set a',
                'results':{'a':3, 'b':4, 'damping':[4,5,6], 'hs':[1,2,3]},
            }
        ]
    compile_results(q)
    compile_damping(q)

if __name__ == '__main__':
    main()