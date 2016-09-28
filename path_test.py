import os

if __name__ == '__main__':
    ABSPATH = os.path.abspath(os.sys.argv[0])
    ABSPATH = os.path.dirname(ABSPATH) + "\\"

    print ABSPATH
    print os.path.abspath('a' + '.xls')
    os.system('pause')
