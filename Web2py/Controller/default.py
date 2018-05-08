# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# This is a sample controller
# this file is released under public domain and you can use without limitations
# -------------------------------------------------------------------------

# ---- example index page ----

def index():

    return dict(message=T('Welcome to VCFWeb!'))

#@auth.requires_login()
def variants():

    import json

    # Select all the records, to show how
    # datatables.net paginates.
    # Rows can't be serialized because they contain a reference to
    # an open database connection. Use as_list()
    # to serialize the query result.
    #data = json.dumps(db(db.TRANSCRIPT).select().as_list())

    data = json.dumps(db.executesql('SELECT consequence, impact, symbol, biotype, sift, polyphen, name, max_freq FROM TRANSCRIPT INNER JOIN PREDICTORS ON TRANSCRIPT.trans_id=PREDICTORS.fk_trans INNER JOIN VARIANT ON VARIANT.var_id=TRANSCRIPT.fk_var INNER JOIN EXISTING_VARIATION ON VARIANT.var_id=EXISTING_VARIATION.fk_var INNER JOIN (select fk_trans,max(freq) as max_freq from (select fk_trans, max(afr_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(amr_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(asj_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(eas_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(fin_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(nfe_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(oth_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(sas_af) as freq from pop_mafs group by fk_trans )a group by fk_trans) b ON b.fk_trans=TRANSCRIPT.trans_id;',as_dict=True))
    #data2 = db.executesql('SELECT name FROM VARIANT INNER JOIN EXISTING_VARIATION ON VARIANT.var_id=EXISTING_VARIATION.fk_var;',as_dict=True)
    #mafs = db.executesql('select max(freq) from (select fk_trans, max(afr_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(amr_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(asj_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(eas_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(fin_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(nfe_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(oth_af) as freq from pop_mafs group by fk_trans union all select fk_trans,  max(sas_af) as freq from pop_mafs group by fk_trans )a;',as_dict=True)
    #data3 = data + str(mafs)
    # Convert to XML for DataTable
    return dict(results=XML(data))


# ---- API (example) -----
@auth.requires_login()
def api_get_user_email():
    if not request.env.request_method == 'GET': raise HTTP(403)
    return response.json({'status':'success', 'email*':auth.user.email})

def hello():
    return dict(message='Hello %(first_name)s' % auth.user)

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

# ---- Action for login/register/etc (required for auth) -----
def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    http://..../[app]/default/user/bulk_register
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    also notice there is http://..../[app]/appadmin/manage/auth to allow administrator to manage users
    """

    return dict(form=auth())


def milogin(): return dict(formulario=auth.login())
def milogout(): return dict(formulario=auth.logout())
def miregistro(): return dict(formulario=auth.register())
def miperfil(): return dict(formulario=auth.profile())

# ---- action to server uploaded static content (required) ---
@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)
