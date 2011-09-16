def run_expression( expr, vars ):
    return eval( expr, {}, vars )

def run_script( expr, vars ):
    exec( expr, {}, vars );
    return locals()['vars']

