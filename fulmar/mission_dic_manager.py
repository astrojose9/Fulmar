import json
from astropy.table import Table

##############################################################################


def write_json_dic(dicname, filename):
    with open(filename, 'w') as output_file:
        output_file.write(json.dumps(dicname))


def read_json_dic(json_dic):
    with open(json_dic) as f:
        t2t = f.read()
    return json.loads(t2t)


def tessDic(tessCat):
    tess_table = Table.read(tessCat, format='ascii.csv', comment='#')
    tess_table.rename_columns(tess_table.colnames, [
                              c.replace(' ', '_')
                              for c in tess_table.colnames])
    tess_table['TIC'] = ['TIC' + x for x in tess_table['tid'].astype(str)]

    tess_table['TOI'] = ['TOI-' + x
                         for x in tess_table['toipfx'].astype(str)]

    tic2toi_dic = dict(zip(tess_table['TIC'], tess_table['TOI']))
    toi2tic_dic = dict(zip(tess_table['TOI'], tess_table['TIC']))

    write_json_dic(tic2toi_dic, 'TIC2TOI.json')
    write_json_dic(toi2tic_dic, 'TOI2TIC.json')

    # with open('TIC2TOI.json', 'w') as tic2toi_file:
    #     tic2toi_file.write(json.dumps(tic2toi_dic))


def keplerDic(keplerCat):
    kep_table = Table.read(keplerCat, format='ascii.csv', comment='#')
    kep_table.rename_columns(kep_table.colnames, [
        c.replace(' ', '_')
        for c in kep_table.colnames])
    kep_table['KIC'] = ['KIC' + x for x in kep_table['kepid'].astype(str)]

    kep_table['kep'] = [x.split(' ')[0]
                        for x in kep_table['kepler_name'].astype(str)]

    KIC2Kepler_dic = dict(zip(kep_table['KIC'], kep_table['kep']))
    Kepler2KIC_dic = dict(zip(kep_table['kep'], kep_table['KIC']))

    write_json_dic(KIC2Kepler_dic, 'KIC2Kepler.json')
    write_json_dic(Kepler2KIC_dic, 'Kepler2KIC.json')

    # with open('KIC2Kepler.json', 'w') as KIC2Kepler_file:
    #     KIC2Kepler_file.write(json.dumps(KIC2Kepler_dic))


def K2Dic(K2Cat):
    K2_table = Table.read(K2Cat, format='ascii.csv', comment='#')
    K2_table.rename_columns(K2_table.colnames, [
        c.replace(' ', '_')
        for c in K2_table.colnames])
    K2_table['EPIC'] = [x.replace(' ', '')
                        for x in K2_table['epic_id'].astype(str)]

    K2_table['K2'] = [x.split(' ')[0]
                      for x in K2_table['k2_name'].astype(str)]

    EPIC2K2_dic = dict(zip(K2_table['EPIC'], K2_table['K2']))
    K22EPIC_dic = dict(zip(K2_table['K2'], K2_table['EPIC']))

    write_json_dic(EPIC2K2_dic, 'EPIC2K2.json')
    write_json_dic(K22EPIC_dic, 'K22EPIC.json')
    # with open('EPIC2Kepler.json', 'w') as EPIC2K2_file:
    #     EPIC2K2_file.write(json.dumps(EPIC2K2_dic))


if __name__ == '__main__':
    tessCat = 'tess_cat.csv'

    keplerCat = 'kepler_cat.csv'

    k2Cat = 'k2_cat.csv'

    tessDic(tessCat)
    keplerDic(keplerCat)
    K2Dic(k2Cat)
