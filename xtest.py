def error(a):
    if a > 10:
        raise ValueError('oops!')
    else:
        print(a)

    return a


error(11)