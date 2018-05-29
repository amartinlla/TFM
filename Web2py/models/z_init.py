# -*- coding: utf-8 -*-
# Description:
# Automatically adds a user to the database called "studio1"

user = 'test2'
pwd = 'test2'
study = 'estudio1'


if not db().select(db.auth_user.ALL).first():
    db.auth_user.insert(
        username = user,
        password = db.auth_user.password.validate(pwd)[0],
        email = 'null@null.com',
        first_name = 'Trial',
        study = study
    )

def hello():
    return dict(message='Hello %(first_name)s' % auth.user)


def goodbye():
    return dict(message='YOU MUST LOG IN FIRST\n')
