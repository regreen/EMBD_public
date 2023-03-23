"""
Function to populate and modify embd table in mySQL database
@author: Becky
"""

import MySQLdb
import json
import numpy as np
from varLookup import PARAM_NAMES, MYSQL_DB_CREDENTIALS


def initialize_db():
    """
    initializes the database
    :return: connection and cursor to the database
    """
	
    myConn = MySQLdb.connect(host=mySQLpath, db=login,
                             read_default_file=MYSQL_DB_CREDENTIALS, read_default_group=login)
    return myConn, myConn.cursor()


def print_rows(cursor, limit=None):
    """
    prints all genes present in the database
    :return:
    """
    sql = "SELECT * FROM genes"
    if limit is not None:
        sql += " LIMIT {}".format(limit)
    cursor.execute(sql)
    rows = cursor.fetchall()
    for row in rows:
        print(row)


def get_columns(cursor, table):
    """
    print table column names
    :return:
    """
    cursor.execute("SHOW COLUMNS FROM {t}".format(t=table))
    result = [value[0] if value[0] != 'bY' else '`bY`' for value in cursor.fetchall()]
    return result


def get_embryo_id(rna, strain, emb_id, date, cursor):
    """
    returns id for fist embryo with given parameters
    :param rna: embd rnai number
    :param strain: strain
    :param emb_id: embryo number (int)
    :param date: date(str)
    :param cursor: mysql cursor
    :return: int
    """
    date = date[0:4] + '-' + date[4:6] + '-' + date[6:8]
    sql = "SELECT id FROM embryos WHERE rnai={r} and strain=\"{s}\" and emb_id={e} and DATE(date)=\"{d}\"".format(r=rna,
                                                                                            s=strain, e=emb_id, d=date)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]


def get_all_ids(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT id FROM genes"
    cursor.execute(sql)
    return list(np.array(cursor.fetchall()).ravel())


def create_param_names(cursor):
    """
    adds parameters to the database
    :return:
    """

    columns = get_columns(cursor, 'genes')
    p_name_prev = ''
    for p_name in PARAM_NAMES:
        print(p_name)
        if p_name + '_GLS' not in columns:
            sql = "ALTER TABLE genes ADD COLUMN {} FLOAT NULL".format(p_name + '_GLS')
            if p_name_prev != '':
                sql += " AFTER {}".format(p_name_prev + '_GLS')
            cursor.execute(sql)
        if p_name + '_MS' not in columns:
            sql = "ALTER TABLE genes ADD COLUMN {} FLOAT NULL".format(p_name + '_MS')
            if p_name_prev != '':
                sql += " AFTER {}".format(p_name_prev + '_MS')
            cursor.execute(sql)


def add_float_column_2embryos(cursor, col):
    """
    adds float column to the embryos database
    :return:
    """
    cursor.execute("SHOW COLUMNS FROM embryos")
    col_names = [value[0] for value in cursor.fetchall()]
    if col not in col_names:
        print('adding {}'.format(col))
        sql = "ALTER TABLE embryos ADD COLUMN {} FLOAT NULL".format(col)
        cursor.execute(sql)


def insert_row(colms, cursor, table):
    """
    adds a new row to specified table in mySQL database
    :param colms: dictionary with column names and values
    :param cursor: cursor
    :param table: db table
    :return: None
    """

    names = ""
    values = ""
    my_colms = get_columns(cursor, table)
    for col in colms:  # cols is the keyname for key value pair in dict
        if col in my_colms and colms[col] is not None:  # colms[col] is the value for key value pair in dict
            if isinstance(colms[col], str):
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=colms[col])
            elif isinstance(colms[col], np.ndarray):
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=json.dumps(colms[col].tolist()))
            elif isinstance(colms[col], (list, dict)):
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=json.dumps(colms[col]))
            elif not np.isnan(colms[col]):  # if a single value (float or int)
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=colms[col])
            elif np.isnan(colms[col]):
                names += '{pn}, '.format(pn=col)
                values += 'NULL, '.format(pv=colms[col])
    names = names[:-2]  # removes the final comma and space from the string
    values = values[:-2]
    sql = "INSERT INTO {t} ({n}) VALUES ({v})".format(t=table, n=names, v=values)
    cursor.execute(sql)


def get_row(myid, cursor, table, colms=None, sonn=False, germ=False, ns=False):
    """
    reads columns in a row from mysql database and returns a dictionary
    :param myid: row id (int)
    :param cursor: mysql cursor
    :param table: table name (str)
    :param colms: list of names of columns
    :return: dictionary
    """

    sql = "SELECT "
    my_colms = get_columns(cursor, table)
    if colms is None:
        sql += '*  '
    else:
        for col in colms:
            if col in my_colms:
                sql += '{c}, '.format(c=col)
    if sonn:
        sql = sql[:-2] + ' FROM {t} WHERE id_sonn_key={id}'.format(t=table, id=myid)
    elif germ:
        sql = sql[:-2] + ' FROM {t} WHERE id_germ_key={id}'.format(t=table, id=myid)
    elif ns:
        sql = sql[:-2] + ' FROM {t} WHERE id_ns_key={id}'.format(t=table, id=myid)
    else:
        sql = sql[:-2] + ' FROM {t} WHERE id={id}'.format(t=table, id=myid)

    cursor.execute(sql)
    vals = list(cursor.fetchall()[0])
    for j in range(len(vals)):
        if isinstance(vals[j], str) and ('[' in vals[j] or '{' in vals[j]):
            vals[j] = json.loads(vals[j])
            if isinstance(vals[j], list):
                vals[j] = np.array(vals[j])
        elif vals[j] is None:
            vals[j] = np.nan
    field_names = [d[0] for d in cursor.description]
    return dict(zip(field_names, vals))


def update_row(myid, colms, cursor, table, sonn=False, germ=False):
    """
    update columns in a given row
    :param myid: row id
    :param colms: dictionary with column names and values
    :param cursor: cursor
    :param table: db table
    :return:
    """

    sql = "UPDATE {t} SET ".format(t=table)
    my_colms = get_columns(cursor, table)
    for col in colms:
        if col in my_colms and colms[col] is not None:
            if isinstance(colms[col], str):
                sql += '{pn}=\'{pv}\', '.format(pn=col, pv=colms[col])
            elif isinstance(colms[col], np.ndarray):
                sql += '{pn}=\'{pv}\', '.format(pn=col, pv=json.dumps(colms[col].tolist()))
            elif isinstance(colms[col], (list, dict)):
                sql += '{pn}=\'{pv}\', '.format(pn=col, pv=json.dumps(colms[col]))
            elif not np.isnan(colms[col]):
                sql += '{pn}={pv}, '.format(pn=col, pv=colms[col])
            elif np.isnan(colms[col]):
                sql += '{pn}=NULL, '.format(pn=col)
    if sonn:
        sql = sql[:-2] + ' WHERE id_sonn_key={id}'.format(id=myid)
    elif germ:
        sql = sql[:-2] + ' WHERE id_germ_key={id}'.format(id=myid)
    else:
        sql = sql[:-2] + ' WHERE id={id}'.format(id=myid)
    cursor.execute(sql)


def add_embryo(rnai, strain, emb_num, date, folder, field, cursor):
    """
    adds embryo to the database
    :param rnai: int
    :param strain: string
    :param emb_num: int
    :param date: string
    :param folder: string
    :param field: int
    :param cursor: mysqldb cursor
    :return: None
    """
    sql = "INSERT INTO embryos (rnai, strain, emb_id, date, file_path, well, field) "
    date = date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' ' + date[9:11] + ':' + date[11:13] + ':' + date[13:15]
    well = int(field[1:3])
    field = int(field[4:])
    sql += "VALUES ({r}, \"{s}\", {e}, \"{d}\", \"{fo}\", {w}, {fi})".format(r=rnai, s=strain, e=emb_num, d=date,
                                                                             fo=folder, w=well, fi=field)
    cursor.execute(sql)


def get_tot_num_embs(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT COUNT(*) FROM embryos"
    cursor.execute(sql)
    return cursor.fetchall()[0][0]


def get_all_params(cursor, p_name, strain, rnai=None):
    """
    returns a set of p_name parameters for all embryos that belong to a given rnai
    :param cursor: mysql cursor
    :param p_name: parameter name
    :param strain: strain name
    :param rnai: rnai value (int). 0 returns values for controls, None returns values for all rnais excluding controls.
    :return: ndarray of parameters
    """
    from emb_handler import get_embryo_revision_version
    version = get_embryo_revision_version()
    if p_name == 'bY':
        p_name = '`bY`'
    if rnai is None:
        sql = "SELECT {pn} FROM embryos WHERE strain=\"{s}\" and dead=0 and version=\"{v}\" and rnai!=0".format(
            pn=p_name, v=version, s=strain)
    else:
        sql = "SELECT {pn} FROM embryos WHERE strain=\"{s}\" and dead=0 and version=\"{v}\" and rnai={r}".format(
            pn=p_name, v=version, s=strain, r=rnai)
    cursor.execute(sql)
    res = []
    for val, in cursor.fetchall():
        if isinstance(val, float) or isinstance(val, int):
            res.append(val)
    return np.array(res)

def get_embryo_number_by_rna(rna, strain='MS'):
    """
    returns list of embryo numbers in mySQL for RNAi condition in MS strain (used to identify and fix head params, tailHead and devHead)
    :param rna: embd rnai number
    :param strain: strain
    :param cursor: mysql cursor
    :return: list of embryo numbers
    """
    conn, cursor = initialize_db()
    sql = "SELECT id, emb_id, rnai FROM embryos WHERE rnai={r} and strain=\"{s}\"".format(r=rna,s=strain)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    conn.close()
    sql_rows_array = np.asarray(sql_rows)[:,1]  # converts tuple to array and returns just the second column, which contains the embryo number
    if len(sql_rows_array) == 0:
        return None
    else:
        return sql_rows_array

def get_file_path_by_rna(rna, strain='MS'):
    """
    returns filepath embryo numbers in mySQL for RNAi condition in MS strain (used to identify and fix head params, tailHead and devHead)
    :param rna: embd rnai number
    :param strain: strain
    :param cursor: mysql cursor
    :return: list of embryo numbers
    """
    conn, cursor = initialize_db()
    sql = "SELECT id, emb_id, rnai, file_path FROM embryos WHERE rnai={r} and strain=\"{s}\"".format(r=rna,s=strain)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    conn.close()
    sql_rows_array = np.asarray(sql_rows)[:,3]  # converts tuple to array and returns just the filepath column, which contains the embryo number
    if len(sql_rows_array) == 0:
        return None
    else:
        return sql_rows_array

def initialize():
    embs = get_embryo_number_by_rna(2, strain='MS')
    print(embs)

if __name__ == '__main__':
    # conn, curs = initialize_db()
    # get_all_params(curs, 'aG', 'MS', 53)
    # conn.close()
    initialize()
