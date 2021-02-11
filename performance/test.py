import argparse
ap = argparse.ArgumentParser("AER302 Performance example")
ap.add_argument('-plane')
ap.add_argument('-V', nargs='?')
ap.add_argument('-H', nargs='?')
ap.add_argument('--endurance', 
				action='store_true',
				default=False,
				help='compute the endurance of the plane')
