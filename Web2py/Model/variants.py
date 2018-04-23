# -*- coding: utf-8 -*-

import psycopg2

db =  DAL('postgres://test2:test2@localhost/estudio1',pool_size=0,lazy_tables=True) # The database is now connected and the connection is stored in the global variable db.

print db._uri # Uniform Resource Identifier
print db._dbname

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
migrate=True)
db.TRANSCRIPT.trans_id.requires = IS_NOT_IN_DB(db, db.TRANSCRIPT.trans_id)
db.TRANSCRIPT.fk_var.requires = IS_NOT_IN_DB(db, db.TRANSCRIPT.fk_var)
db.TRANSCRIPT.fk_var.requires = IS_IN_DB(db, db.VARIANT.var_id, '%(nombre)s')

db.define_table('PREDICTORS',
Field('fk_trans','integer'),
Field('sift'),
Field('polyphen'),
Field('loftool'),
Field('cadd_phred'),
Field('cadd_raw'),
migrate=True)
db.PREDICTORS.fk_trans.requires = IS_NOT_IN_DB(db, db.PREDICTORS.fk_trans)
db.PREDICTORS.fk_trans.requires = IS_IN_DB(db, db.TRANSCRIPT.trans_id, '%(nombre)s')



#for reg in db().select(db.VARIANT.ALL):
#    print reg.var_id

#print db.executesql('SELECT * FROM VARIANT;')

# If no database exists, generate a database of 101 unique records
# with names in the form John1 Smith1, John43 Smith43, etc.
#if db(db.variant).isempty():
#    print 'variant table is empty'
