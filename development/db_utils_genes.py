"""
Function to populate and modify embd table in mySQL database
"""

import pandas as pd
import db_utils_embryos
from fileLookup import FILENAME_MANUAL_SCORING, FILENAME_GENE_WB_INFO, FILENAME_EMBD_INFO, FILENAME_GENE_CANCER_RANK
from embdFunc import getEMBDNumber



def get_columns(cursor, table):
    """
    print table column names
    :return:
    """
    cursor.execute("SHOW COLUMNS FROM {t}".format(t=table))
    result = [value[0] for value in cursor.fetchall()]
    return result


def print_rows_ce_id(id_list, cursor):
    """
    prints all genes present in the database with ce_id from list
    :param id_list: list if genes sequences
    :return:
    """
    if isinstance(id_list, str):
        id_list = [id_list]
    sql = "SELECT * FROM genes WHERE ce_id in ("
    for ce_id in id_list:
        sql += "\"{}\",".format(ce_id)
    sql = sql[:-1] + ")"
    cursor.execute(sql)
    rows = cursor.fetchall()
    for row in rows:
        print row


def check_gene_present(gene_name, cursor):
    """
    checks if the gene with given id present in the database
    :param gene_name: gene id
    :return: bool
    """
    cursor.execute("SELECT * FROM GENES WHERE gene_name='{}'".format(gene_name))
    rows = cursor.fetchall()
    return len(rows) > 0


def get_cancer_rank(h_gene, cr_info):
    """
    populates cancer rank information from csv file containing human gene names and cancer rank based on Divoli et al
    :param takes in human gene name (i.e. GATA3)
    :return: the highest ranked value for the gene between the tumor suppressor or oncogene cancer ranking list
    """
    
    ts_ind, og_ind = cr_info.loc[cr_info['Tumor Sup'] == h_gene], cr_info.loc[cr_info['Oncogene'] == h_gene]
#     ts_ind, og_ind = np.where(cr_info[:,1]==h_gene), np.where(cr_info[:,2]==h_gene)
    if not ts_ind.empty or not og_ind.empty:
        rank = min(ts_ind.iloc[0,0], og_ind.iloc[0,0])
        return rank
    else:
        return None

    
def update_cancer_rank(cursor):
    """
    updates cancer rank information- looks up human gene names and finds and applies the minimum rank value (i.e. highest cancer rank) for each list of orthologs
    :param 
    :return:
    """
    
    cr_info = read_cancer_rank()
    for my_id in db_utils_embryos.get_all_ids(curs):
        sql = "SELECT id, h_gene_name, cancer_rank FROM genes WHERE id='{}'".format(my_id)
        cursor.execute(sql)
        id, h_gene_name, cancer_rank = cursor.fetchall()[0]
        if h_gene_name is not None:
            h_gene_name_list = h_gene_name.split(',')
            rank_list = []
            for h_gene in h_gene_name_list:
                h_rank = get_cancer_rank(h_gene, cr_info)
                if h_rank is not None:
                    rank_list.append(h_rank)
            if len(rank_list)>0:
                rank = min(rank_list)
                sql = "UPDATE genes SET cancer_rank=\"{cr}\" WHERE id={id}".format(cr=rank, id=my_id)
                cursor.execute(sql)

def update_human_ortho2_germ(embd, curs):
    '''gets human ortholog- looks up ortholog from Ortholist 2, and updates an existing row in mySQL'''
    from fileLookup import FILENAME_EMBD_ORTHOLIST
    import pandas as pd

    sql = "SELECT id, embd_id, ce_id, h_gene_name FROM genes WHERE embd_id={}".format(embd)
    curs.execute(sql)
    id, embd_id, ce_id, h_gene_name = curs.fetchall()[0]

    file = FILENAME_EMBD_ORTHOLIST
    data = pd.read_csv(file, ',')
    self_data = data[data['Locus ID'] == ce_id]
    self_data_sort = self_data.sort_values('No. of Databases', ascending=False)
    h_gene_list = self_data_sort['HGNC Symbol'].values.tolist()
    h_gene_str = ", ".join([str(elem) for elem in h_gene_list])
    print(h_gene_str)
    h_genes = h_gene_str

    sql2 = "UPDATE genes SET h_gene_name=\"{h}\" WHERE id={id};".format(h=h_genes, id=id)
    # sql2 = "UPDATE sonnichsen SET top_CSI_matches=\'{h}\' WHERE id_sonn_key={id};".format(h=json.dumps(hit_list), id=self.id_sonn_key)

    curs.execute(sql2)
    conn.commit()


def update_embd_info():
    """
    populates or updates mySQL genes table with information from csv file containing EMBD number, gene names,
    oligo/plate info, emb lethal data, etc
    :param 
    :return:
    """
 
    # embd_info = read_embd_info()
    # for index, row in embd_info.iterrows():
    #     if row['gene target'] != None: #and not check_gene_present(row['Public Name'])
    #         sql = "SELECT id, embd_id, oligo_1, oligo_2, oligo_coord, oligo_plate, date_imaged, imaged_by," \
    #               "imaging_coord, room temp, start_temp, end_temp, EMB_GLS, EMB_MS, abL1_GLS, abL1_MS, FROM genes WHERE ce_id='{sn}'".format(sn=row['gene target'])
    #         cursor.execute(sql)
    #         sql_rows = cursor.fetchall()
    #         embd_id_present = False
    #         for s_row in sql_rows:
    #             if s_row[1] is None:
    #                 sql = "UPDATE genes SET "
    #                 if row['gene target'] != 'N.A.':
    #                     sql += "embd_id=\"{}\",".format(int(row['embd_number'][-4:]))
    #                 if row['T7 oligo'] != 'N.A.':
    #                     sql += "oligo_1=\"{}\",".format(row['T7 oligo'])
    #                 if row['T3 oligo'] != 'N.A.':
    #                     sql += "oligo_2=\"{}\",".format(row['T3 oligo'])
    #                 if row['Coordinate'] != 'N.A.':
    #                     sql += "oligo_coord=\"{}\",".format(row['Coordinate'])
    #                 if row['plate name'] != 'N.A.':
    #                     sql += "oligo_plate=\"{}\",".format(row['plate name'])
    #                 if row['Date Imaged'] != 'N.A.':
    #                     sql += "date_imaged=\"{}\",".format(row['Date Imaged'])
    #                 if row['imaged_by'] != 'N.A.':
    #                     sql += "imaged_by=\"{}\",".format(row['imaged_by'])
    #                 if row['imaging_coord'] != 'N.A.':
    #                     sql += "imaging_coord=\"{}\",".format(row['imaging_coord'])
    #                 sql = sql[:-1] + " WHERE id={id}".format(id=s_row[0]) #sql[:-1] removes the comma from the last line
    #                 cursor.execute(sql)
    #                 embd_id_present = True
    #             elif s_row[0] == row['embd_number']:
    #                 embd_id_present = True

    conn, cursor = db_utils_embryos.initialize_db()
    table = "genes"
    embd_info = read_embd_info()
    for index, row in embd_info.iterrows():
        if row['gene target'] != None:
            print("updating {0}").format(row['embd_number'])
            id = get_id_from_embd(row['embd_number'], cursor)  # this is the mySQL unique identifier for the given EMBD number
            columns = {}  # creates a dictionary for each id to be passed to update_row function
            if row['gene target'] != 'N.A.':
                columns['ce_id'] = row['gene target']
            if row['T7 oligo'] != 'N.A.':
                columns['oligo_1'] = row['T7 oligo']
            if row['T3 oligo'] != 'N.A.':
                columns['oligo_2'] = row['T3 oligo']
            if row['Coordinate'] != 'N.A.':
                columns['oligo_coord'] = row['Coordinate']
            if row['plate name'] != 'N.A.':
                columns['oligo_plate'] = row['plate name']
            if row['RNA_concentration'] != 'N.A.':
                columns['RNA_concentration'] = row['RNA_concentration']
            if row['Date_Imaged'] != 'N.A.':
                columns['date_imaged'] = row['Date_Imaged']
            if row['imaged_by'] != 'N.A.':
                columns['imaged_by'] = row['imaged_by']
            if row['imaging_coord'] != 'N.A.':
                columns['imaging_coord'] = row['imaging_coord']
            if row['room_temp'] != 'N.A.':
                columns['room_temp'] = row['room_temp']
            if row['start_temp'] != 'N.A.':
                columns['start_temp'] = row['start_temp']
            if row['end_temp'] != 'N.A.':
                columns['end_temp'] = row['end_temp']
            if row['EMB_GLS'] != 'N.A.':
                columns['EMB_GLS'] = row['EMB_GLS']
            if row['EMB_MS'] != 'N.A.':
                columns['EMB_MS'] = row['EMB_MS']
            if row['abL1_GLS'] != 'N.A.':
                columns['abL1_GLS'] = row['abL1_GLS']
            if row['abL1_MS'] != 'N.A.':
                columns['abL1_MS'] = row['abL1_MS']
            db_utils_embryos.update_row(id, columns, cursor, table)
    conn.commit()
    conn.close()

def update_interacting_genes():
    '''updates the interacting genes for each gene in the database by reading from wb_info csv file.
    Interacting genes were identified by wormbase queries'''
    conn, cursor = db_utils_embryos.initialize_db()
    table = "genes"
    wb_info = read_wb_info()
    for index, row in wb_info.iterrows():
        if row['Your Input'] != None:
            gene = row['Your Input']
            embd_number = getEMBDNumber(gene)
            print("updating {0}").format(embd_number)
            id = get_id_from_genename(gene, cursor)  # this is the mySQL unique identifier for the given gene
            columns = {}  # creates a dictionary for each id to be passed to update_row function
            if row['Interacting_Gene'] != 'N.A.':
                columns['wb_interacting_genes'] = row['Interacting_Gene']
                db_utils_embryos.update_row(id, columns, cursor, table)
    conn.commit()
    conn.close()

def read_manual_data():
    """
    reads information from csv file that contains all manually scored data
    :return: info
    """
    manual_data = pd.read_csv(FILENAME_MANUAL_SCORING)
    return manual_data

def update_manual_data(cursor):
    """
    populates genes with information from csv file containing manual scoring data by EMBD number. To update all entries,
     'replace' must be set to true. Otherwise, only new entries (with no populated data) will be changed.
    :param
    :return:
    """
    replace = False # change to True to re-populate all data. Leave as false to update only new empty rows
    manual_data = read_manual_data()
    for index, row in manual_data.iterrows():
        if row['embd_id'] != None:
            sql = "SELECT id, embd_id, man_GSembs_dors, man_GSembs_vent, man_GSembs_lat, man_GSembs_tot, man_MSembs_dors," \
                  " man_MSembs_vent, man_MSembs_lat, man_MSembs_tot, man_GS1, man_GS2, man_GS3, man_GS4, man_GS5, man_GS6," \
                  " man_GS7, man_GS8, man_GS9, man_GS12, man_GS14, man_GS15, man_GS19, man_GS20, man_GS21, " \
                  "man_GS22, man_GSemb_comments, man_MS1, man_MS2, man_MS4, man_MS5, man_MS6, man_MS7, man_MS8, " \
                  "man_MS9, man_MS10, man_MS11, man_MS13, man_MS14, man_MS15, man_MS17, man_MS19, man_MS20, man_MS21," \
                  " man_MS22, man_MSemb_comments FROM genes WHERE embd_id='{sn}'".format(sn=row['embd_id'])
            cursor.execute(sql)
            sql_rows = cursor.fetchall()
            manual_info_present = False
            for s_row in sql_rows:
                if s_row[2] is None or replace:
                # if s_row[2] is None or s_row[1]<7:
                    print(s_row[0])  # prints SQL unique ID
                    sql = "UPDATE genes SET "
                    if row['man_GSembs_dors'] != 'N.A.':
                        sql += "man_GSembs_dors=\"{}\",".format(int(row['man_GSembs_dors']))
                    if row['man_GSembs_vent'] != 'N.A.':
                        sql += "man_GSembs_vent=\"{}\",".format(int(row['man_GSembs_vent']))
                    if row['man_GSembs_lat'] != 'N.A.':
                        sql += "man_GSembs_lat=\"{}\",".format(int(row['man_GSembs_lat']))
                    if row['man_GSembs_tot'] != 'N.A.':
                        sql += "man_GSembs_tot=\"{}\",".format(int(row['man_GSembs_tot']))
                    if row['man_MSembs_dors'] != 'N.A.':
                        sql += "man_MSembs_dors=\"{}\",".format(int(row['man_MSembs_dors']))
                    if row['man_MSembs_vent'] != 'N.A.':
                        sql += "man_MSembs_vent=\"{}\",".format(int(row['man_MSembs_vent']))
                    if row['man_MSembs_lat'] != 'N.A.':
                        sql += "man_MSembs_lat=\"{}\",".format(int(row['man_MSembs_lat']))
                    if row['man_MSembs_tot'] != 'N.A.':
                        sql += "man_MSembs_tot=\"{}\",".format(int(row['man_MSembs_tot']))
                    if row['man_GS1'] != 'N.A.':
                        sql += "man_GS1=\"{}\",".format(int(row['man_GS1']))
                    if row['man_GS2'] != 'N.A.':
                        sql += "man_GS2=\"{}\",".format(int(row['man_GS2']))
                    if row['man_GS3'] != 'N.A.':
                        sql += "man_GS3=\"{}\",".format(int(row['man_GS3']))
                    if row['man_GS4'] != 'N.A.':
                        sql += "man_GS4=\"{}\",".format(int(row['man_GS4']))
                    if row['man_GS5'] != 'N.A.':
                        sql += "man_GS5=\"{}\",".format(int(row['man_GS5']))
                    if row['man_GS6'] != 'N.A.':
                        sql += "man_GS6=\"{}\",".format(int(row['man_GS6']))
                    if row['man_GS7'] != 'N.A.':
                        sql += "man_GS7=\"{}\",".format(int(row['man_GS7']))
                    if row['man_GS8'] != 'N.A.':
                        sql += "man_GS8=\"{}\",".format(int(row['man_GS8']))
                    if row['man_GS9'] != 'N.A.':
                        sql += "man_GS9=\"{}\",".format(int(row['man_GS9']))
                    if row['man_GS12'] != 'N.A.':
                        sql += "man_GS12=\"{}\",".format(int(row['man_GS12']))
                    if row['man_GS14'] != 'N.A.':
                        sql += "man_GS14=\"{}\",".format(int(row['man_GS14']))
                    if row['man_GS15'] != 'N.A.':
                        sql += "man_GS15=\"{}\",".format(int(row['man_GS15']))
                    if row['man_GS19'] != 'N.A.':
                        sql += "man_GS19=\"{}\",".format(int(row['man_GS19']))
                    if row['man_GS20'] != 'N.A.':
                        sql += "man_GS20=\"{}\",".format(int(row['man_GS20']))
                    if row['man_GS21'] != 'N.A.':
                        sql += "man_GS21=\"{}\",".format(int(row['man_GS21']))
                    if row['man_GS22'] != 'N.A.':
                        sql += "man_GS22=\"{}\",".format(int(row['man_GS22']))
                    if row['man_MS1'] != 'N.A.':
                        sql += "man_MS1=\"{}\",".format(int(row['man_MS1']))
                    if row['man_MS2'] != 'N.A.':
                        sql += "man_MS2=\"{}\",".format(int(row['man_MS2']))
                    if row['man_MS10'] != 'N.A.':
                        sql += "man_MS10=\"{}\",".format(int(row['man_MS10']))
                    if row['man_MS4'] != 'N.A.':
                        sql += "man_MS4=\"{}\",".format(int(row['man_MS4']))
                    if row['man_MS5'] != 'N.A.':
                        sql += "man_MS5=\"{}\",".format(int(row['man_MS5']))
                    if row['man_MS6'] != 'N.A.':
                        sql += "man_MS6=\"{}\",".format(int(row['man_MS6']))
                    if row['man_MS7'] != 'N.A.':
                        sql += "man_MS7=\"{}\",".format(int(row['man_MS7']))
                    if row['man_MS8'] != 'N.A.':
                        sql += "man_MS8=\"{}\",".format(int(row['man_MS8']))
                    if row['man_MS9'] != 'N.A.':
                        sql += "man_MS9=\"{}\",".format(int(row['man_MS9']))
                    if row['man_MS11'] != 'N.A.':
                        sql += "man_MS11=\"{}\",".format(int(row['man_MS11']))
                    if row['man_MS13'] != 'N.A.':
                        sql += "man_MS13=\"{}\",".format(int(row['man_MS13']))
                    if row['man_MS14'] != 'N.A.':
                        sql += "man_MS14=\"{}\",".format(int(row['man_MS14']))
                    if row['man_MS15'] != 'N.A.':
                        sql += "man_MS15=\"{}\",".format(int(row['man_MS15']))
                    if row['man_MS17'] != 'N.A.':
                        sql += "man_MS17=\"{}\",".format(int(row['man_MS17']))
                    if row['man_MS19'] != 'N.A.':
                        sql += "man_MS19=\"{}\",".format(int(row['man_MS19']))
                    if row['man_MS20'] != 'N.A.':
                        sql += "man_MS20=\"{}\",".format(int(row['man_MS20']))
                    if row['man_MS21'] != 'N.A.':
                        sql += "man_MS21=\"{}\",".format(int(row['man_MS21']))
                    if row['man_MS22'] != 'N.A.':
                        sql += "man_MS22=\"{}\",".format(int(row['man_MS22']))
                    if row['man_GSemb_comments'] != 'N.A.':
                        sql += "man_GSemb_comments=\"{}\",".format(row['man_GSemb_comments'])
                    if row['man_MSemb_comments'] != 'N.A.':
                        sql += "man_MSemb_comments=\"{}\",".format(row['man_MSemb_comments'])
                    sql = sql[:-1] + " WHERE id={id}".format(id=s_row[0])  # sql[:-1] removes the comma from the last line
                    print(sql)
                    cursor.execute(sql)
                    manual_info_present = True

                elif s_row[2] == row['man_GSembs_dors']:
                    manual_info_present = True
                    print('manual info present already')
                else:
                    print('wrong column')



def update_gene_wb_info(cursor):
    """
    populates genes with information from wormbase
    :param id_list: list if genes to update their ce names
    :return:
    """

    wb_info = read_wb_info()
    for index, row in wb_info.iterrows():
        if row['Status'] == 'Live':
            sql = "SELECT id, gene_name FROM genes WHERE ce_id='{sn}'".format(sn=row['Your Input'])
            cursor.execute(sql)
            sql_rows = cursor.fetchall()
            gene_present = False
            for s_row in sql_rows:
                if s_row[1] is None:
                    sql = "UPDATE genes SET "
                    if row['Public Name'] != 'N.A.':
                        sql += "gene_name=\"{}\",".format(row['Public Name'])
                    if row['WormBase Gene ID'] != 'N.A.':
                        sql += "wbgene_id=\"{}\",".format(row['WormBase Gene ID'])
                    if row['Description Text'] != 'N.A.':
                        if row['Description Type']=='Automated_description':
                            sql += "ce_description=\"(A) {}\",".format(row['Description Text'])
                        else:
                            sql += "ce_description=\"(C) {}\",".format(row['Description Text'])
                    if row['Human Ortholog'] != 'N.A.':
                        sql += "ensembl_id=\"{}\",".format(row['Human Ortholog'])
                    if row['Disease Info'] != 'N.A.':
                        sql += "disease_info=\"{}\",".format(row['Disease Info'])
                    if row['Uniprot'] != 'N.A.':
                        sql += "uniprot=\"{}\",".format(row['Uniprot'])
                    if row['TreeFam'] != 'N.A.':
                        sql += "treefam=\"{}\",".format(row['TreeFam'])
                    if row['RNAi Phenotype Observed'] != 'N.A.':
                        sql += "rnai_ph_obs=\"{}\",".format(row['RNAi Phenotype Observed'])
                    if row['RNAi Phenotype Not Observed'] != 'N.A.':
                        sql += "rnai_ph_notobs=\"{}\",".format(row['RNAi Phenotype Not Observed'])
                    if row['Allele Phenotype Observed'] != 'N.A.':
                        sql += "allele_ph_obs=\"{}\",".format(row['Allele Phenotype Observed'])
                    if row['Allele Phenotype Not Observed'] != 'N.A.':
                        sql += "allele_ph_notobs=\"{}\",".format(row['Allele Phenotype Not Observed'])
                    if row['Expr_pattern Tissue'] != 'N.A.':
                        sql += "expr_tissue=\"{}\",".format(row['Expr_pattern Tissue'])
                    if row['Expr_pattern LifeStage'] != 'N.A.':
                        sql += "expr_stage=\"{}\",".format(row['Expr_pattern LifeStage'])
                    if row['Interacting_Gene'] != 'N.A.':
                        sql += "wb_interacting_genes=\"{}\",".format(row['Interacting_Gene'])
                    sql = sql[:-1] + " WHERE id={id}".format(id=s_row[0])
                    cursor.execute(sql)
                    gene_present = True
                elif s_row[1] == row['Public Name']:
                    gene_present = True
            if not gene_present:
                sql = "INSERT INTO genes "
                col_str = "(ce_id,"
                val_str = " VALUES (\"{}\",".format(row['Your Input'])
                if row['Public Name'] != 'N.A.':
                    col_str += "gene_name,"
                    val_str += "\"{}\",".format(row['Public Name'])
                if row['WormBase Gene ID'] != 'N.A.':
                    col_str += "wbgene_id,"
                    val_str += "\"{}\",".format(row['WormBase Gene ID'])
                if row['Description Text'] != 'N.A.':
                    if row['Description Type']=='Automated_description':
                        col_str += "ce_description,"
                        val_str += "\"(A) {}\",".format(row['Description Text'])
                    else:
                        col_str += "ce_description,"
                        val_str += "\"(C) {}\",".format(row['Description Text'])
                if row['Human Ortholog'] != 'N.A.':
                    col_str += "ensembl_id,"
                    val_str += "\"{}\",".format(row['Human Ortholog'])
                if row['Disease Info'] != 'N.A.':
                    col_str += "disease_info,"
                    val_str += "\"{}\",".format(row['Disease Info'])
                if row['Uniprot'] != 'N.A.':
                    col_str += "uniprot,"
                    val_str += "\"{}\",".format(row['Uniprot'])
                if row['TreeFam'] != 'N.A.':
                    col_str += "treefam,"
                    val_str += "\"{}\",".format(row['TreeFam'])
                if row['RNAi Phenotype Observed'] != 'N.A.':
                    col_str += "rnai_ph_obs,"
                    val_str += "\"{}\",".format(row['RNAi Phenotype Observed'])
                if row['RNAi Phenotype Not Observed'] != 'N.A.':
                    col_str += "rnai_ph_notobs,"
                    val_str += "\"{}\",".format(row['RNAi Phenotype Not Observed'])
                if row['Allele Phenotype Observed'] != 'N.A.':
                    col_str += "allele_ph_obs,"
                    val_str += "\"{}\",".format(row['Allele Phenotype Observed'])
                if row['Allele Phenotype Not Observed'] != 'N.A.':
                    col_str += "allele_ph_notobs,"
                    val_str += "\"{}\",".format(row['Allele Phenotype Not Observed'])
                if row['Expr_pattern Tissue'] != 'N.A.':
                    col_str += "expr_tissue,"
                    val_str += "\"{}\",".format(row['Expr_pattern Tissue'])
                if row['Expr_pattern LifeStage'] != 'N.A.':
                    col_str += "expr_stage,"
                    val_str += "\"{}\",".format(row['Expr_pattern LifeStage'])
                if row['Interacting_Gene'] != 'N.A.':
                    col_str += "wb_interacting_genes,"
                    val_str += "\"{}\",".format(row['Interacting_Gene'])
                sql = sql + col_str[:-1] + ") " + val_str[:-1] + ")"
                cursor.execute(sql)


# def update_ens_refs(cursor): FIXME
#     """
#     updates the gene names from ensembl ids
#     Note: this will not run on windows due to issue with pyensembl dependence on unix only resource package
#     :return:
#     """
#
#     for my_id in get_all_ids():
#         sql = "SELECT ensembl_id, h_gene_name FROM genes WHERE id='{}'".format(my_id)
#         cursor.execute(sql)
#         ens_string, h_n_assigned = cursor.fetchall()[0]
#         if ens_string is not None:
#             all_ens_ids = [ens_id.split('|')[0] for ens_id in ens_string.split(',')]
#             h_names = get_ensemble_ref(all_ens_ids)
#             all_h_genes = ','.join(h_names)
#             if all_h_genes != h_n_assigned:
#                 sql = "UPDATE genes SET h_gene_name=\"{h}\" WHERE id={id}".format(h=all_h_genes, id=my_id)
#                 cursor.execute(sql)


def get_unppl_id(col):
    pass


def remove_gene(id, cursor):
    cursor.execute("DELETE FROM genes WHERE embd_id={}".format(id))


def read_wb_info():
    """
    reads information provided by the wormbase.org
    :return: info
    """
    wb_info = pd.read_csv(FILENAME_GENE_WB_INFO)
    return wb_info


def read_embd_info():
    """
    reads information provided by csv- contains the EMBD assignment, gene, oligos, oligo plate # and position
    :return: info
    """
    embd_info = pd.read_csv(FILENAME_EMBD_INFO)
    return embd_info


def read_cancer_rank():
    """
    reads information provided by csv- contains human gene names and their cancer ranks based on Divoli et al
    :return: info
    """
    cr_info = pd.read_csv(FILENAME_GENE_CANCER_RANK)
    return cr_info


def get_tot_num_genes(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT COUNT(*) FROM genes"
    cursor.execute(sql)
    return cursor.fetchall()[0][0]


def get_id_from_embd(embd_id, cursor):
    """
    returns id for first gene with embd number
    :param embd_id: embd number
    :param cursor: db cursor
    :return: int unique identifier for mySQL table
    """

    if isinstance(embd_id, basestring):
        embd_id = int(embd_id[-4:])
    sql = "SELECT id FROM genes WHERE embd_id={}".format(embd_id)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]

def get_id_from_sonn(sonn_id, curs):
    sql = "SELECT id_sonn_key FROM sonnichsen WHERE sonn_id={}".format(sonn_id)
    curs.execute(sql)
    sql_rows = curs.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]

def get_id_from_germ(germ_id, curs):
    sql = "SELECT id_germ_key FROM germline WHERE germ_id={}".format(germ_id)
    curs.execute(sql)
    sql_rows = curs.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]

def get_id_from_ns(ns_id, curs):
    sql = "SELECT id_ns_key FROM nuclear_screen WHERE ns_id={}".format(ns_id)
    curs.execute(sql)
    sql_rows = curs.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]

def get_id_from_genename(genename, cursor):
    """
    returns id for first gene with genename (i.e.'R05H5.2')
    :param embd_id: genename
    :param cursor: db cursor
    :return: int unique identifier for mySQL table
    """

    sql = "SELECT id FROM genes WHERE ce_id=\"{}\"".format(genename)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]


def get_max_id(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT id FROM genes ORDER BY id DESC LIMIT 0, 1"
    cursor.execute(sql)
    return cursor.fetchall()[0][0]


def get_ensemble_ref(ens_name_list):
    import pyensembl
    ensembl = pyensembl.EnsemblRelease()
    if not isinstance(ens_name_list, list):
        ens_name_list = [ens_name_list]
    names = [ensembl.gene_by_id(ens_name).gene_name for ens_name in ens_name_list]
    return names



def create_columns_pad_rank(cursor):
    """
    adds columns to pad_rank database for top 50 hits (rank#) and pad values for the top 50 hits (pad#)
    - use this to build the pad_rank database
    :return:
    """

    columns = get_columns(cursor, 'pad_rank')  # retrieves columns from table as list
    for i in range(1, 51):  # all potential RNAi conditions in the EMBD screen
        print(i)
        if i not in columns:
            sql = "ALTER TABLE pad_rank ADD COLUMN `{}` INT(11) NULL".format('rank' + '{}'.format(i))  # adds columns for the top 50 RNAi
            # conditions (result formatted as INT)
            cursor.execute(sql)
            sql = "ALTER TABLE pad_rank ADD COLUMN `{}` FLOAT NULL".format('pad' + '{}'.format(i))# adds columns for PAD values for the
            # top 50 RNAi conditions (result formatted as float)
            cursor.execute(sql)

def get_embryo_count_by_RNAi(embd_id, cursor):
    """
    returns number of embryos in SQL for GLS and MS by embd number
    :param embd_id: embd number
    :param cursor: db cursor
    :return: #GLS embs, [list of GLS embs], #MS embs, [list of MS embs]
    """

    if isinstance(embd_id, basestring):
        embd_id = int(embd_id[-4:])
    sql = "SELECT count(emb_id) FROM embryos WHERE strain LIKE 'GLS' AND rnai={}".format(embd_id);
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        GLS_count = None
    else:
        GLS_count = sql_rows[0][0]
    sql = "SELECT count(emb_id) FROM embryos WHERE strain LIKE 'MS' AND rnai={}".format(embd_id);
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        MS_count = None
    else:
        MS_count = sql_rows[0][0]
    # print(int(GLS_count),int(MS_count))
    return(GLS_count, MS_count)


def get_reciprocal_matches(seed, top_hit_number, cursor, thr=0.55):
    """
    Gets top PAD hits from SQL and checks to see if the hit is reciprocal. Criteria is top x hits in first direction, top
    x hits in reverse direction... as long as PAD values are above 0.55 (or whatever thresh is established)
    :param seed: EMBD number
    :param top_hit_number: the number of top hits to consider in first direction (x)
    :param thresh- adjustable to determine what threshold of PAD to consider for best hits
    :return: best hit list (a list of RNAi conditions that show reciprocal top hits)
    """
    import numpy as np

    ranks_considered = []  # generates a list of 'rank1, rank2, rank3...' based on the specified top hit number
    pad_considered = []  # generates a list of 'pad1, pad2, pad3...' based on top hit number

    for i in range(1, top_hit_number+1):  # generates the two lists above
        rank = 'rank{}'.format(i)
        ranks_considered.append(rank)
        pad_val = 'pad{}'.format(i)
        pad_considered.append(pad_val)

    ranks_considered_string = ','.join(ranks_considered)  # converts ranks_considered to a string for SQL query
    pad_considered_string = ','.join(pad_considered)  # converts pad_considered to a string for SQL query


    sql = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(ranks_considered_string, seed)
    sql2 = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(pad_considered_string, seed)

    cursor.execute(sql)
    hits = cursor.fetchall()
    hits = np.array(hits[0])  # convert retrieved dataframe to np.array
    cursor.execute(sql2)
    pad = cursor.fetchall()
    pad = np.array(pad[0])

    high_pad = pad > thr  # generates a TRUE/FALSE list to ID usable indices
    inds = np.where(high_pad==True)  # retrieves inds where high_pad is true to apply to ranks considered
    ranks_use = hits[inds]  # a new list of ranks to use in the query that have PADs >0.5

    recip_hit_number = top_hit_number  # generates a reciprocal rank value (if top hit number is 10, this will be also be 10)
    # recip_hit_number = top_hit_number * 2  # generates a reciprocal rank value (if top hit number is 10, this will be 20)
    recip_ranks_considered = []  # generates a list of 'rank1, rank2, rank3...' based on the specified top hit number

    for j in range(1, recip_hit_number+1):  # generates the two lists above #This was missing the +1 in the first iteration- is fixed
        recip_rank = 'rank{}'.format(j)
        recip_ranks_considered.append(recip_rank)
    recip_ranks_considered_string = ','.join(recip_ranks_considered)

    best_hits = []  # list of hits that are reciprocal hits. Populated by loop below
    for hit in ranks_use:
        sql3 = "SELECT {0} FROM PAD_RANK WHERE rnai={1};".format(recip_ranks_considered_string, hit)
        cursor.execute(sql3)
        recip_hits = cursor.fetchall()
        recip_hits = recip_hits[0]
        if seed in recip_hits:
            best_hits.append(hit)
    print((seed, best_hits))
    return(best_hits)

def update_automated_description():
    '''
    Populates automated description column based on proximity of gene of interest to centroid of manual groups
    :param cursor:
    :return:
    '''
    from groupAssesment import get_all_closest_groups

    get_all_closest_groups(to_csv=False, to_SQL=True, to_plot=False)


def update_top_hits_PAD(cursor):
    '''
    Populates top_PAD_hits column in mySQL with top 10 PAD hits. All hits must have a PAD above 0.5.
    :param cursor:
    :return:
    '''
    import numpy as np
    from embdFunc import getGeneCommonName
    thr = 0.5

    # for i in range(1,504):
    for i in [1406,674,952,866,885,703,1360,1904,1903,1905]:

        print('starting {0}').format(i)
        seed = i
        ranks_considered = 'rank1, rank2, rank3, rank4, rank5, rank6, rank7, rank8, rank9, rank10'
        pad_considered = 'pad1, pad2, pad3, pad4, pad5, pad6, pad7, pad8, pad9, pad10'
        sql = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(ranks_considered, seed)
        sql2 = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(pad_considered, seed)

        cursor.execute(sql)  # gets EMBD top 10
        hits = cursor.fetchall()
        hits = np.array(hits[0])  # convert retrieved dataframe to np.array

        cursor.execute(sql2)  # gets PADs, top 10
        pad = cursor.fetchall()
        pad = np.array(pad[0])

        high_pad = pad > thr  # generates a TRUE/FALSE list to ID usable indices (above 0.5 PAD)
        inds = np.where(high_pad==True)  # retrieves inds where high_pad is true to apply to ranks considered
        ranks_use = hits[inds]  # a new list of ranks to use in the query that have PADs >thresh

        genes = []  # populated by common gene name for each EMBD_id
        for item in ranks_use:
            genes.append(getGeneCommonName(item))
        genes_use = ', '.join(genes)  # makes a string to insert into mySQL query

        sql = "SELECT id, embd_id, top_hits_PAD FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, top_hits_PAD = cursor.fetchall()[0]
        sql2 = "UPDATE genes SET top_hits_PAD=\"{h}\" WHERE id={id};".format(h=genes_use, id=id)
        cursor.execute(sql2)
        conn.commit()


def update_WT_like_or_phenotype(cursor):
    '''
    Populates column in SQL to define if RNAi condition is called as 'WT-like' (moving and sustaining movement), or
    'phenotype' (movement is 0 for either GLS, MS or both). Reads movementGLS and movementMS to determine.
    :param cursor:
    :return: populates mySQL
    '''
    for i in range(1,504):
        print('starting {0}').format(i)
        sql = "SELECT id, embd_id, movement_GLS, movement_MS, auto_description, WT_like_or_phenotype FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, movementG, movementM, auto_description, WT_like_or_phenotype = cursor.fetchall()[0]
        if movementG == 0 or movementM == 0:
            sql2 = "UPDATE genes SET WT_like_or_phenotype=\"{p}\" WHERE id={id};".format(p='Phenotype', id=id)
            cursor.execute(sql2)
            conn.commit()
        elif movementG == 1 and movementM == 1:
            if "Embryos develop at a slow pace." in auto_description:
                sql2 = "UPDATE genes SET WT_like_or_phenotype=\"{p}\" WHERE id={id};".format(p='Dev. Delay', id=id)
                cursor.execute(sql2)
                conn.commit()
            else:
                sql2 = "UPDATE genes SET WT_like_or_phenotype=\"{p}\" WHERE id={id};".format(p='WT-like', id=id)
                cursor.execute(sql2)
                conn.commit()


def update_column_reciprocal_matches(cursor):
    """
    populates reciprocal matches column in genes database
    :return:
    """
    from embdFunc import getGeneCommonName

    for i in range(310,311):
        best_hits = get_reciprocal_matches(i, 10, cursor, thr=0.55)
        top_five = best_hits[0:10]  # changed to top 10

        genes = []
        for item in best_hits:
            genes.append(getGeneCommonName(item))
        top_five_genes = ', '.join(genes[0:10])
        print((i,top_five, top_five_genes))

        sql = "SELECT id, embd_id, movement_GLS, movement_MS, WT_like_or_phenotype, reciprocal_matches FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, movementG, movementM, WT_like_or_phenotype, reciprocal_matches = cursor.fetchall()[0]

        if movementG == 0 or movementM == 0:  # checks to see if there is a scored phenotype by checking movement
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=top_five_genes, id=id)
            cursor.execute(sql2)
            conn.commit()
        elif movementG == 1 and movementM == 1 and "Dev. Delay" in WT_like_or_phenotype:  # looks for dev. delay
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=top_five_genes, id=id)
            cursor.execute(sql2)
            conn.commit()
        else:  # if no scored phenotype (movement) or dev delay, doesnt provide 'hits' as it is expected to be wild-type.
            wt_note = 'Embryos do not exhibit an arrest phenotype'
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=wt_note, id=id)
            cursor.execute(sql2)
            conn.commit()

def update_other_names(cursor): # in progress
    '''Reads Other Gene Names from csv file where other names have been pulled from wormbase and populates them into 'other_names' column in mySQL'''

    from fileLookup import FILENAME_OTHER_NAMES, SKIPPED
    import pandas as pd
    count=0
    for i in range(1,1900):
        count+=1
        if 'EMBD{0:04}'.format(i) not in SKIPPED:
            embd = i
            sql = "SELECT id, embd_id, ce_id, other_name FROM genes WHERE embd_id={}".format(embd) #selects key and column from MySQL
            curs.execute(sql)
            id, embd_id, ce_id, other_name = curs.fetchall()[0]
            file = FILENAME_OTHER_NAMES
            data = pd.read_csv(file, ',') # loads other names data from a csv (pulled from Simple Mine)
            row_data = data.loc[data['EMBD_ID']== embd_id, ['gene_key', 'EMBD_ID','Other Name']] #
            if row_data['gene_key'].iloc[0] == id:
                sql2 = "UPDATE genes SET other_name=\"{k}\" WHERE id={id};".format(k=row_data['Other Name'].iloc[0], id=id)
            curs.execute(sql2)
            conn.commit()
            if count%10==0:
                print(count)



def get_sql_data_for_man_group_conditions():
    '''
    retrieves the automated parameters and manual parameters for the conditions selected in manual groupings from the SQL. To be used for classifier algorithms.
    :param cursor:
    :return: df of automated and manual params for 94 manually grouped conditions
    '''
    import numpy as np
    conn, cursor = db_utils_embryos.initialize_db()
    sql = "SELECT id, embd_id, man_group, movement_GLS, movement_MS, WT_like_or_phenotype, maxG_GLS, aG_GLS, sY_GLS, fracR_GLS, fracG_GLS, CoM0Rtd_GLS_GLS, CoM1Rscale_GLS_GLS, MoI0GavgRes_GLS_GLS, MoI1Gscale_GLS_GLS, " \
          "MoI1GstdTail_GLS_GLS, MoI0Rtd_GLS_GLS, MoI0RavgRes_GLS_GLS, MoI1Rtd_GLS_GLS, MoI0Ytd_GLS_GLS, MoI0YavgRes_GLS_GLS, MoI1Yscale_GLS_GLS, MoI1YavgRes_GLS_GLS, sG_MS, maxR_MS, aR_MS, sR_MS, maxHead_MS, tailHead_MS, " \
          "devHead_MS, devLength_MS, tailLength_MS, CoM0Gtd_MS_MS, CoM0Gscale_MS_MS, CoM1Gscale_MS_MS, CoM1GavgRes_MS_MS, CoM0Rscale_MS_MS, CoM0RavgRes_MS_MS, CoM1Rtd_MS_MS, CoM1RavgRes_MS_MS, CoM1RstdTail_MS_MS, " \
          "MoI0Gscale_MS_MS, MoI1Gscale_MS_MS, MoI1GstdTail_MS_MS, MoI1Rscale_MS_MS, MoI1RavgRes_MS_MS, MoI1RstdTail_MS_MS,man_GS1, man_GS2, man_GS3, man_GS4, man_GS5, man_GS6,man_GS7, man_GS8, man_GS9, man_GS12, man_GS14," \
          "man_GS15, man_GS19, man_GS20, man_GS21, man_GS22, man_GSemb_comments, man_MS1, man_MS2, man_MS4, man_MS5, man_MS6, man_MS7, man_MS8, man_MS9, man_MS10, man_MS11, man_MS13, man_MS14, man_MS15, man_MS17, man_MS19," \
          "man_MS20, man_MS21, man_MS22, man_MSemb_comments FROM genes WHERE man_group IS NOT NULL"
    cursor.execute(sql)
    df = pd.read_sql_query(sql, conn)  #note that fetchall doesnt retrieve a df and this command works better
    return df

def force_version(cursor, pv):
    """
    forces version change in mySQL
    :param cursor: mysqldb cursor
    :param pv: version to reset
    :return: None"""
    from emb_handler import get_embryo_revision_version
    from UT_bench_mark_emb import TestBenchCase
    import unittest
    # test = unittest.TestLoader().loadTestsFromTestCase(TestBenchCase)
    # test_result = unittest.TestResult()
    # test.run(test_result)
    # if test_result.wasSuccessful():
    sql = "UPDATE embryos SET version = \"{v}\" WHERE version = \"{pv}\"".format(v=get_embryo_revision_version(), pv=pv)
    cursor.execute(sql)
    print("versions are successfully updated")
    # else:
    #     print("Didn't update version because bench mark test failed")
    #     raise ValueError


if __name__ == '__main__':
    print("starting")
    conn, curs = db_utils_embryos.initialize_db()

    # update_manual_data(curs)
#     test_list = [5, 7, 9, 10, 11, 15, 25, 27, 28, 30, 31, 32, 35, 45, 51, 55, 59, 61, 69, 71, 75, 79, 84, 85, 87, 90, 94, 95, 96, 97, 99, 108, 109, 111, 118, 122, 126, 130, 131, 135, 137, 147, 152, 155,167, 168, 177, 178, 179, 180, 182, 198, 207, 210, 225, 233, 235, 240, 241, 251, 255, 264, 267, 274, 275, 277, 280, 282, 288, 289, 291, 294, 299, 312, 315, 325, 329, 330, 331, 332, 333, 336, 359, 360, 361, 366, 368, 369, 372, 378, 379, 383, 384, 386, 392, 394, 397, 398, 407, 410, 411, 412, 414, 422, 435, 438, 440, 441, 450, 454, 455, 456, 457, 458, 463, 478, 479, 481, 483, 487, 488, 492, 496, 497, 499, 503]
# #this is a list of genes that are in the 10th rank position- used to flag reciprocal match issue
#     fix_list_recip_match = [27, 31, 35, 51, 55, 59, 61, 69, 75, 85, 90, 94, 95, 96, 97, 99, 108, 109, 111, 122, 126, 130, 131, 135, 137, 147, 155, 167, 168, 178, 179, 180, 182, 207, 210, 225, 233, 235, 240, 241, 255, 264, 274, 275, 282, 288, 289, 294, 299, 312, 315, 325, 329, 330, 331, 332, 333, 336, 359, 360, 361, 366, 368, 372, 379, 384, 386, 392, 394, 397, 407, 410, 411, 435, 441, 450, 454, 456, 457, 478, 479, 488, 496, 499, 503]
# these are the flagged conditions that need to be updated on mySQL

    '''get reciprocal hits'''
    # for gene in test_list:
    # genelist = [339,497,35,429,476,201,310]
    # for gene in genelist:
    #     get_reciprocal_matches(gene, 10, curs, thr=0.55)  # args are EMBD#, # of top hits considered, curs,
    # # and pad threshold to use
    # conn.close()

    '''update reciprocal hits'''
    # update_column_reciprocal_matches(curs)
    # conn.close()

    '''populate WT-like column'''
    # update_WT_like_or_phenotype(curs)
    # conn.close()
    #
    '''populate top_PAD_hits column'''
    # update_top_hits_PAD(curs)
    # conn.close()

    '''populate auto_description column SQL'''
    # update_automated_description()

    '''update ortholist2 genes for genes table'''
    # for i in range(1,2000):
    #     embd = i
    #     update_human_ortho2_germ(embd, curs)

#     print(get_embryo_id(1, 'GLS', 1, curs))
    # update_gene_wb_info()
    # update_ens_refs()

    # update_interacting_genes()
    # update_embd_info()
    # print("done")

    # update_cancer_rank(curs)

    '''update other names'''
    update_other_names(curs)
    conn.close()


    # create_param_names()
    # print(get_columns(curs))
    # print_rows(curs, 10)
    # print('starting manual update')
    # update_manual_data(curs)
    # create_columns_pad_rank(curs)
    # add_genes_pad_rank(curs)
    # print(get_id_from_embd(496, curs))
    # print(get_gene_id('EMBD0496', curs))
    # GLS_count, MS_count = get_embryo_count_by_RNAi(3, curs)

    # force_version(curs, 'f6668a2')
    # df = get_sql_data_for_man_group_conditions()
    # conn.commit()
    # conn.close()
