def test_sessionfinish(session, exitstatus):
    """ whole test run finishes. """
    for filename in glob.glob('atesa_v2/tests/test_temp/*'):
        print('\n')
        print(sys.path[0] + '/' + filename)
        os.remove(sys.path[0] + '/' + filename)
