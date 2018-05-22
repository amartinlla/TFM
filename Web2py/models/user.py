from gluon.tools import Auth


db.define_table('USERS',
Field('user_id','integer'),
Field('username'),
Field('password','password'),
Field('study'),
migrate=True)
db.USERS.user_id.requires = IS_NOT_IN_DB(db, db.USERS.user_id)
