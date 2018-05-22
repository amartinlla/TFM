# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# This is a sample controller
# this file is released under public domain and you can use without limitations
# -------------------------------------------------------------------------
import logging
from gluon.storage import Storage


logger = logging.getLogger("web2py.app.vcfweb")
logger.setLevel(logging.DEBUG)


def index():
    #if auth.user.username == session.auth.user.username:
    if auth.is_logged_in():
        redirect(URL('variants'))
        #redirect(URL('..','VCFWeb/filters'))
    else:
        return dict(message=T('Welcome to web2py!'))

def user():
    return dict(form=auth())

def next():
    return dict()


@auth.requires_login()
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


@auth.requires_login()
def filters():
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

# ---- API (example) -----
@auth.requires_login()
def api_get_user_email():
    if not request.env.request_method == 'GET': raise HTTP(403)
    return response.json({'status':'success', 'email*':auth.user.email})


# ---- Smart Grid (example) -----
@auth.requires_membership('admin') # can only be accessed by members of admin groupd
def grid():
    response.view = 'generic.html' # use a generic view
    tablename = request.args(0)
    if not tablename in db.tables: raise HTTP(403)
    grid = SQLFORM.smartgrid(db[tablename], args=[tablename], deletable=False, editable=False)
    return dict(grid=grid)

# ---- Embedded wiki (example) ----
def wiki():
    auth.wikimenu() # add the wiki to the menu
    return auth.wiki() 



# ---- action to server uploaded static content (required) ---
@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)
