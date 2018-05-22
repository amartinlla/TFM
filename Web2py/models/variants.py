# -*- coding: utf-8 -*-
# try something like
#@auth.requires_login()
def index():
    import json

    # Select all the records, to show how
    # datatables.net paginates.
    # Rows can't be serialized because they contain a reference to
    # an open database connection. Use as_list()
    # to serialize the query result.
    #data = json.dumps(db(db.TRANSCRIPT).select().as_list())
    data = json.dumps(db.executesql('SELECT * FROM AUX2 INNER JOIN POP_MAFS ON POP_MAFS.fk_trans=AUX2.trans_id;',as_dict=True))

    # Convert to XML for DataTable
    return dict(results=XML(data))

def secret():
   conditions = [request.client == '127.0.0.1', auth.user]
   if not any(conditions):
       redirect(URL('user', args='login'))
   return dict(message='Hello test2 again!')

#@auth.requires_login()
def variants():

    import json

    # Select all the records, to show how
    # datatables.net paginates.
    # Rows can't be serialized because they contain a reference to
    # an open database connection. Use as_list()
    # to serialize the query result.
    #data = json.dumps(db(db.TRANSCRIPT).select().as_list())
    data = json.dumps(db.executesql('SELECT * FROM AUX2 INNER JOIN POP_MAFS ON POP_MAFS.fk_trans=AUX2.trans_id;',as_dict=True))

    # Convert to XML for DataTable
    return dict(results=XML(data))
