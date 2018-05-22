# -*- coding: utf-8 -*-

import psycopg2


db =  DAL('postgres://test2:test2@localhost/estudio1',pool_size=0,lazy_tables=True) # The database is now connected and the connection is stored in the global variable db.

print db._uri # Uniform Resource Identifier
print db._dbname
from gluon.tools import Auth, Crud, Service, PluginManager


auth = Auth(db)
crud, service, plugins = Crud(db), Service(), PluginManager()

## create all tables needed by auth if not custom tables
auth.define_tables(username=True)

## configure auth policy
auth.settings.actions_disabled.append('register')
auth.settings.registration_requires_verification = True
auth.settings.registration_requires_approval = True
auth.settings.reset_password_requires_verification = True


db.define_table('VARIANT',
Field('var_id','integer'),
Field('chrom','integer'),
Field('pos','integer'),
Field('id'),
Field('ref'),
Field('alt'),
Field('qual'),
Field('filter'),
migrate=True)
db.VARIANT.var_id.requires = IS_NOT_IN_DB(db, db.VARIANT.var_id)
db.VARIANT.chrom.requires = IS_NOT_IN_DB(db, db.VARIANT.chrom)
db.VARIANT.pos.requires = IS_NOT_IN_DB(db, db.VARIANT.pos)

db.define_table('EXISTING_VARIATION',
Field('id','integer'),
Field('fk_var','integer'),
Field('name'),
migrate=True)
db.EXISTING_VARIATION.id.requires = IS_NOT_IN_DB(db, db.EXISTING_VARIATION.id)
db.EXISTING_VARIATION.fk_var.requires = IS_NOT_IN_DB(db, db.EXISTING_VARIATION.fk_var)
db.EXISTING_VARIATION.fk_var.requires = IS_IN_DB(db, db.VARIANT.var_id, '%(nombre)s')

db.define_table('TRANSCRIPT',
Field('trans_id','integer'),
Field('fk_var','integer'),
Field('allele'),
Field('consequence'),
Field('impact'),
Field('symbol'),
Field('gene'),
Field('feature_type'),
Field('feature'),
Field('biotype'),
Field('CLIN_SIG'),
migrate=True)
db.TRANSCRIPT.trans_id.requires = IS_NOT_IN_DB(db, db.TRANSCRIPT.trans_id)
db.TRANSCRIPT.fk_var.requires = IS_NOT_IN_DB(db, db.TRANSCRIPT.fk_var)
db.TRANSCRIPT.fk_var.requires = IS_IN_DB(db, db.VARIANT.var_id, '%(nombre)s')

db.define_table('PREDICTORS',
Field('pred_id','integer'),
Field('fk_var','integer'),
Field('fk_trans','integer'),
Field('sift'),
Field('polyphen'),
Field('loftool'),
Field('cadd_phred'),
Field('cadd_raw'),
migrate=True)
db.PREDICTORS.pred_id.requires = IS_NOT_IN_DB(db, db.PREDICTORS.pred_id)
db.PREDICTORS.fk_var.requires = IS_NOT_IN_DB(db, db.PREDICTORS.fk_var)
db.PREDICTORS.fk_var.requires = IS_IN_DB(db, db.VARIANT.var_id, '%(nombre)s')
db.PREDICTORS.fk_trans.requires = IS_NOT_IN_DB(db, db.PREDICTORS.fk_trans)
db.PREDICTORS.fk_trans.requires = IS_IN_DB(db, db.TRANSCRIPT.trans_id, '%(nombre)s')


db.define_table('POP_MAFS',
Field('mafs_id','integer'),
Field('fk_var','integer'),
Field('fk_trans','integer'),
Field('max_af'),
Field('af'),
Field('afr_af'),
Field('amr_af'),
Field('asj_af'),
Field('eas_af'),
Field('fin_af'),
Field('nfe_af'),
Field('oth_af'),
Field('sas_af'),
migrate=True)
db.POP_MAFS.mafs_id.requires = IS_NOT_IN_DB(db, db.POP_MAFS.mafs_id)
db.POP_MAFS.fk_var.requires = IS_NOT_IN_DB(db, db.POP_MAFS.fk_var)
db.POP_MAFS.fk_var.requires = IS_IN_DB(db, db.VARIANT.var_id, '%(nombre)s')
db.POP_MAFS.fk_trans.requires = IS_NOT_IN_DB(db, db.POP_MAFS.fk_trans)
db.POP_MAFS.fk_trans.requires = IS_IN_DB(db, db.TRANSCRIPT.trans_id, '%(nombre)s')

'''
db.define_table('OTHER_INFO',
Field('other_id','integer'),
Field('fk_trans','integer'),
Field('exon'),
Field('intron'),
Field('HGVSc'),
Field('HGVSp'),
Field('cDNA'),
Field('CDS'),
Field('proton'),
Field('aminoacids'),
Field('codons'),
Field('distance'),
Field('strand'),
Field('flags'),
Field('variant_class'),
Field('symbol_source'),
Field('HGNC_ID'),
Field('CCDS'),
Field('ENSP'),
Field('swissprot'),
Field('trembl'),
Field('uniparc'),
Field('gene_pheno'),
Field('domains'),
Field('HGVS_OFFSET'),
Field('SOMATIC'),
Field('PHENO'),
Field('PUBMED'),
migrate=True)
db.OTHER_INFO.other_id.requires = IS_NOT_IN_DB(db, db.OTHER_INFO.other_id)
db.OTHER_INFO.fk_trans.requires = IS_NOT_IN_DB(db, db.OTHER_INFO.fk_trans)
db.OTHER_INFO.fk_trans.requires = IS_IN_DB(db, db.TRANSCRIPT.trans_id, '%(nombre)s')
'''

db.define_table('AUX2',
Field('var_id','integer'),
Field('trans_id','integer'),
Field('chrom','integer'),
Field('pos','integer'),
Field('consequence'),
Field('impact'),
Field('symbol'),
Field('gene'),
Field('feature_type'),
Field('feature'),
Field('biotype'),
Field('CLIN_SIG'),
Field('name'),
Field('sift'),
Field('polyphen'),
Field('loftool'),
Field('cadd_phred'),
migrate=True)
db.AUX2.var_id.requires = IS_NOT_IN_DB(db, db.AUX2.var_id)
db.AUX2.trans_id.requires = IS_NOT_IN_DB(db, db.AUX2.trans_id)

db.define_table('USERS',
Field('user_id','integer'),
Field('username'),
Field('password','password'),
Field('study'),
migrate=True)
db.USERS.user_id.requires = IS_NOT_IN_DB(db, db.USERS.user_id)



#for reg in db().select(db.VARIANT.ALL):
#    print reg.var_id

#print db.executesql('SELECT * FROM VARIANT;')

# If no database exists, generate a database of 101 unique records
# with names in the form John1 Smith1, John43 Smith43, etc.
#if db(db.variant).isempty():
#    print 'variant table is empty'
