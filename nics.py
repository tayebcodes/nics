'''
AUTHOR: TAYEB KAKESHPOUR
'''
import numpy, math, argparse, ast, sys, os
import pybel, rdkit
from rdkit import Chem
from pprint import pprint
from scipy.optimize import curve_fit
from subprocess import Popen, PIPE, call

'''
How to run this script:

make sure there is Gaussian input file names "3c.com" exists in the directore and tpye:

python nics.py -I 3c.com -R "1 2 3 4 5 6"

'''
# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input", help="input file name: it can be eithe Gaussian16 input (*.com) or output (.log)")
parser.add_argument("-R", "--rings", help="ring element indices starting from 1. Also, makes sure single quotes are included. e.g: '1 3 6 7 9'")
parser.add_argument("-M", "--method", help="method: e.g mpwpw91")
parser.add_argument("-B", "--basis", help="basis set: e.g 6-311+G(2d,p)")
parser.add_argument("-C", "--charge", help="charge : e.g -1")
parser.add_argument("-S", "--spin", help="basis set: e.g 2")
args = parser.parse_args()

input = args.input
rings = args.rings
method = args.method
basis = args.basis
charge = args.charge
spin = args.spin

# set default value for level of theory if not specified.
if method == None:
	method = 'mpwpw91'

if basis == None:
	basis = '6-311+(2d,p)'
if charge == None:
	charge = '0'
if spin == None:
	spin = '1'

gaussian_card = '# nmr=giao gen nosymm ' + method
charge_and_multiplicity = charge + ' ' + spin

atomlist = {'H':1,'He':2,\
			'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,\
			'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,\
			'K':19,'Ca':20,'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,\
			'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,\
			'Rb':37,'Sr':38,'In':49,'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,\
			'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,\
			'Cs':55,'Ba':56,'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,\
			'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,\
			'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,\
			'Fr':87,'Ra':88,\
			'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,\
			'Rf':104,'Db':15,'Sg':106,'Bh':107,'Hs':108,'Mt':109,}

###########################################################################################################
########################################## DEFINING FUNCTIONS #############################################
###########################################################################################################
# functions for getting geomtries:
def getGeometryGeneral(input):
	'''
	this function tries to distingush between a gaussian input and output and extract geometries from them.
	'''
	if input.endswith('.log'):
		result = get_lines(input,'Standard orientation','Rotational constants (GHZ)',5,-1,'1 3 4 5')
	elif input.endswith('.com'):
		result = getGINgeometry(input)
	elif input.endswith('.mol2'):
		result = None
	elif input.endswith('.xyz'):
		result = None
	elif input.endswith('.mol'):
		result = None
	elif input.endswith('.sdf'):
		result = None
	elif input.endswith('.pdb'):
		result = None
	elif input.endswith('.xml'):
		result = None
	return result

def getGINgeometry(input):
	'''
	This function extracts the xyz coordinates from a xyz format Gaussian16 input file.
	it might work with Gaussian09 and Gaussian03 input files, too. At this point, the 
	extention of the file does not matter.
	'''
	# open the file read it in and close it.
	with open(input) as f:
		gin = f.readlines()
	# split each line into list to make a nested list.
	lines = [i.split() for i in gin]
	# convert the numbers in the list to floats if possible.
	final = []
	for i in range(len(lines)):
		final0 = []
		for j in range(len(lines[i])):
			b = lines[i][j]
			try:
				a = float(b)
			except ValueError:
				# print("error")
				a = b
			# print(type(a))
			final0.append(a)
		final.append(final0)
	# check the lines tha have atomic coordinates and add them to the final list.
	lines = final
	geometry = []
	for line in lines:
		if len(line) == 4 and type(line[1]) == float and type(line[2]) == float and type(line[3]) == float:
			geometry.append(line)
	return geometry

def get_lines(g16,start_phrase,end_phrase,lines_after_start,lines_before_end,columns_matrix):
	'''
	This function gets the line between two phrases and extracts the requested columns between them.
	Also, lines before and after the end and start phrase can be eliminated.
	start_phrase is where the point of interest begins.
	lines_after_atart is the number of lines to be eliminated after the start phrase.
	end phrase is the point where marks the point after which we are not interested in.
	lines_before_end is the number of lines that need to be eliminated before the end phrase
	columns_matrix is a string of integers separated by space which show the columns of interest starting from 0.
	'''
	# Open the file and read it as a list by each line
	with open(g16) as f:
	    gout = f.readlines()
	# Find the lines that have the geometries in between
	start_phrases = []
	for i, elem in enumerate(gout):
	    if start_phrase in elem:
	        start_phrases.append(i)
	end_phrases = []
	for i, elem in enumerate(gout):
	    if end_phrase in elem:
	        end_phrases.append(i)
	# Get the last geometry in the file.
	start_line = start_phrases[-1] + lines_after_start
	end_line = end_phrases[-2] + lines_before_end
	# Sometime the frequencies cause confusion so correct for that.
	if (start_line > end_line):
		end_line = end_phrases[-1] + lines_before_end
	# Split each atom's properties and delete useless columns from the 2D list.
	lines = gout[start_line:end_line]
	lines = [i.split() for i in lines]
	columns_matrix0 = columns_matrix.split()
	columns = [int(i) for i in columns_matrix0]
	rows = list(range(0,(len(lines))))
	final = []
	for row in rows:
		final0 = []
		for column in columns:
			b = lines[row][column]
			a = ast.literal_eval(b)
			final0.append(a)
		final.append(final0)
	#pprint(final)
	return final

# functions for finding mean plane:
def SortByAtomicNumber(unsorted):
	'''
	This function sort a list of atoms by their atomic numbers.
	'''
	for i, atom in enumerate(unsorted):
		a = unsorted[i][0]
		b = atomlist.get(a)
		atom.insert(0,b)
	sorted1 = sorted(unsorted, key=lambda x: x[0], reverse=True)

	for i, atom in enumerate(unsorted):
		atom.pop(0)
	# print2D(sorted1,False)
	return sorted1

def getAxis(geometry,column):
	'''
	this function extracts a column of interest from a 2D nested list.
	'''
	if column == 'X' or column == 'x':
		column = 1
	elif column == 'Y' or column == 'y':
		column = 2
	elif column == 'Z' or column == 'z':
		column = 3
	axis = []
	for i, row in enumerate(geometry):
		element = geometry[i][column]
		axis.append(element)
		# print(axis)
	return axis

def getRings(geometry):
	'''
	gets the ring elements matrix using rdkit. both openbabel and rdkit need to be installed. 
	'''
	# write an xyz file of the original coordinates
	orig_stdout = sys.stdout
	number = len(geometry)
	sys.stdout = open('temp.xyz','wt')
	print(number)
	print('This is a comment line.')
	print2D(geometry,False)
	sys.stdout.close()
	sys.stdout=orig_stdout
	# convert the xyz file into a mol file to be opened by rdkit
	call(['obabel', 'temp.xyz', '-O', 'temp.mol'], stdout=PIPE, stderr=PIPE)
	# Popen(['obabel', '1.xyz', '-O', '1.mol'], stdout=PIPE, stderr=PIPE)	
	m = Chem.MolFromMolFile('temp.mol')
	ssr = Chem.GetSymmSSSR(m)
	# call(['rm', 'temp.xyz', 'temp.mol'])
	return ssr

def getRingAxis(centered_geometry,ring_elements,column):
	'''
	this function extract partial columns of rinterest from a 2D nested list of geometry.
	'''
	if column == 'X' or column == 'x':
		column = 1
	elif column == 'Y' or column == 'y':
		column = 2
	elif column == 'Z' or column == 'z':
		column = 3
	ring_elements = ring_elements.split()
	ring_elements = [int(i) for i in ring_elements]
	ring_elements = [i - 1 for i in ring_elements]
	# print(ring_elements)
	axis = getAxis(centered_geometry,column)

	ring_axis = []

	for row in ring_elements:
		element = centered_geometry[row][column]
		ring_axis.append(element)
	return ring_axis

def centerRing(geometry,ring_elements,x_column,y_column,z_column):
	ring_elements = ring_elements.split()
	ring_elements = [int(i) for i in ring_elements]
	ring_elements = [i - 1 for i in ring_elements]
	# print(ring_elements)
	x = getAxis(geometry,x_column)
	y = getAxis(geometry,y_column)
	z = getAxis(geometry,z_column)
	x_ring = []
	y_ring = []
	z_ring = []
	for row in ring_elements:
		element = geometry[row][1]
		x_ring.append(element)
	for row in ring_elements:
		element = geometry[row][2]
		y_ring.append(element)
	for row in ring_elements:
		element = geometry[row][3]
		z_ring.append(element)	

	x_center = sum(x_ring)/len(x_ring)
	y_center = sum(y_ring)/len(y_ring)
	z_center = sum(z_ring)/len(z_ring)
	x_centered = [round(element - x_center,10) for element in x]
	y_centered = [round(element - y_center,10) for element in y]
	z_centered = [round(element - z_center,10) for element in z]
	centered_geometry = []
	for i,element in enumerate(x_centered):
		atom = []
		e = geometry[i][0]
		atom.append(e)
		x = x_centered[i]
		atom.append(x)
		y = y_centered[i]
		atom.append(y)
		z = z_centered[i]
		atom.append(z)
		centered_geometry.append(atom)
	return centered_geometry

def findMeanPlane(geometry,ring_elements):
	'''
	Given a 2D nested list of xyz ccordinates (geomtery) and a string containing the indices of the ring elements
	in th form of 'int int int int ...' (ring_elements, a string separated by spaces) finds the a, b, and c of a plane.
	The geometry does not need to be centered
	'''
	# center the geometry and get the ring atom coordinates
	centered_geometry = centerRing(geometry,ring_elements,'x','y','z')
	x_ring = getRingAxis(centered_geometry,ring_elements,'x')
	y_ring = getRingAxis(centered_geometry,ring_elements,'y')
	z_ring = getRingAxis(centered_geometry,ring_elements,'z')
	# find the best fit
	def function(big_x,a,b):
		'''
		this is a function for ring plane.
		'''
		x,y = big_x
		return (a*x) + (b*y)

	p0 = 1, 1
	best_vals, covar = curve_fit(function, (x_ring,y_ring), z_ring, p0)
	[a0,b0,c0] = [best_vals[0], best_vals[1], -1]
	# convert the mean plane vector to a unit vector
	length = (a0**2 + b0**2 + c0**2)**(0.5)
	[a,b,c] = [-a0/length, -b0/length, -c0/length]
	return (a,b,c)

def insertNICS(geometry,ring_elements):
	'''
	This function takes a 2D nested list of atomic coordinates, center them to the ring_elemens which is 
	a string of integers separated by spaces, finds the mean plane, and rotates the molecule so the vector oof
	the mean plane is in thhe Z direction. It also add both NICS(1)zz and NICS(0) probes
	'''
	atoms = getAxis(geometry,0)  # get atoms column
	centered_geometry = centerRing(geometry,ring_elements,'x','y','z') # center the ring 
	[a,b,c] = findMeanPlane(centered_geometry,ring_elements)
	[rho,tetha,phi] = car2sph(a,b,c) # convert the mean plane from cartesian to spherical coordinates

	xx = getAxis(centered_geometry,'x')
	yy = getAxis(centered_geometry,'y')
	zz = getAxis(centered_geometry,'z')
	centered = [atoms,xx,yy,zz]
	#print2D(centered,True)
	##################################### ROTATE AROUND Z BY -THETA ######################################
	[xx1,yy1,zz1] = TethaRot(xx,yy,zz,tetha) # rotate the molecule arounz Z by -tetha
	[aa1,bb1,cc1] = TethaRot(a,b,c,tetha) # rotate the plane vector arounz Z by -tetha
	[rho1,tetha1,phi1] = car2sph(aa1,bb1,cc1) # convert again the plane vector into spherical coordinates
	# print('rotated by tetha around z rho =>  tetha phi = ', rho1, tetha1, phi1)
	# tetha_around_z = [atoms,xx1,yy1,zz1]
	# print2D(tetha_around_z, True)
	##################################### ROTATE AROUND Y BY -PHI ########################################
	[xx2,yy2,zz2] = PhiRot(xx1,yy1,zz1,-phi1) # rotate the molecule around Y by phi
	[aa2,bb2,cc2] = PhiRot(aa1,bb1,cc1,-phi1) # rotate the plane vector around Y by phi
	# print('rotated by tetha around y by phi => rho tetha phi = ', rho1, tetha1, phi1)
	# phi_around_y = [atoms,xx2,yy2,zz2]
	# print2D(phi_around_y, True)
	##################################### ROTATE AROUND Y BY PHI ########################################
	[xx3,yy3,zz3] = PhiRot(xx1,yy1,zz1,phi1) # rotate the molecule around Y by phi
	[aa3,bb3,cc3] = PhiRot(aa1,bb1,cc1,phi1) # rotate the plane vector around Y by phi
	# print('rotated by tetha around y by -phi => rho tetha phi = ', rho1, tetha1, phi1)
	# minus_phi_around_y = [atoms,xx3,yy3,zz3]
	# print2D(minus_phi_around_y, True)
	# find which courdinates to use
	if abs(cc2) > abs(cc3):
		[xx_final,yy_final,zz_final] = [xx2,yy2,zz2]
	else:
		[xx_final,yy_final,zz_final] = [xx3,yy3,zz3]
	xx_final  = [ '%.8f' % elem for elem in xx_final ]
	yy_final  = [ '%.8f' % elem for elem in yy_final ]
	zz_final  = [ '%.8f' % elem for elem in zz_final ]
	xx_final = numpy.array(xx_final).tolist()
	yy_final = numpy.array(yy_final).tolist()
	zz_final = numpy.array(zz_final).tolist()
	# format the columns and get them ready for printing.
	xx_final = [i.rjust(14,' ') for i in xx_final]
	yy_final = [i.rjust(14,' ') for i in yy_final]
	zz_final = [i.rjust(14,' ') for i in zz_final]
	# add the nics probes and transpose the matrix
	nics_bq = [['bq','    0.00000000','    0.00000000','    0.00000000'],\
	['bq','    0.00000000','    0.00000000','    1.00000000'],\
	['bq','    0.00000000','    0.00000000','   -1.00000000']]
	nics = [atoms,xx_final,yy_final,zz_final]
	nics = [list(i) for i in zip(*nics)]
	nics.append(nics_bq[0])
	nics.append(nics_bq[1])
	nics.append(nics_bq[2])
	return nics

# functions for manipulating the coordinates:
def car2sph(x,y,z):
	'''
	This function converts cartestion coordinates of a xyz point to spherical cooridinates.
	'''
	rho   = (x**2+y**2+z**2)**(0.5)
	tetha = -math.atan(y/x)
	phi   = -math.atan(((x**2+y**2)**0.5)/z)
	return (rho, tetha, phi)

def TethaRot(x,y,z,tetha):
	'''
	This function rotates coordinates of a single point or a set of points around the Z-AXIS by theta.
	'''
	x0 = numpy.array(x)
	y0 = numpy.array(y)
	z0 = numpy.array(z)
	x1 = (x0*math.cos(tetha))-(y0*math.sin(tetha))
	y1 = (x0*math.sin(tetha))+(y0*math.cos(tetha))
	z1 = z0
	return (x1,y1,z1)

def PhiRot(x,y,z,phi):
	'''
	This function rotates coordinates of a single point or a set of points around the Y-AXIS by phi.
	'''
	x0 = numpy.array(x)
	y0 = numpy.array(y)
	z0 = numpy.array(z)
	x1 = (x0*math.cos(phi))+(z0*math.sin(phi))
	y1 =  y0
	z1 = (z0*math.cos(phi))-(x0*math.sin(phi))
	return (x1,y1,z1)

# printing tool:
def print2D(list_name,transpose):
	''' This fuction prints a 2D list in a nice way with columns and everything lined up
	the list name is the name of the 2D list that you want to print. the transpose part is for wether
	you want to transpose the mateix before printing. if so type "True", otherwise type "False".
	'''
	if transpose == True:
		list_name = [list(i) for i in zip(*list_name)]
	mx = max((len(str(ele)) for sub in list_name for ele in sub))
	for row in list_name:
		print(" ".join(["{:<{mx}}".format(ele,mx=mx) for ele in row]))

def delete_file(file_name):
    """
    Deletes the specified file if it exists in the current directory.
    Handles the case where the file doesn't exist gracefully.
    
    Args:
        file_name (str): The name of the file to delete.
    """
    try:
        # Try to delete the file
        os.remove(file_name)
    except FileNotFoundError:
        # Handle the case where the file doesn't exist gracefully
        pass
    except Exception as e:
        # Handle any other exceptions that might occur
        pass


###########################################################################################################
########################################## END OF FUNCTIONS ###############################################
###########################################################################################################

# import geometry from a Gaussian input or output.
unsorted_geometry = getGeometryGeneral(input)

# check if ring atoms are requested manually.
if rings == None:
	# if ring atoms are not specified, sort the atoms and find the rings.
	geometry = SortByAtomicNumber(unsorted_geometry)
	ssr = getRings(geometry)
else:
	# otherwise, do not sort the geometry just convert the ring indices that were input to an appropriate nested list.
	geometry = unsorted_geometry
	ssr = rings.split()
	ssr = [int(x) for x in ssr]
	ssr = [x -1 for x in ssr]
	ssr = [ssr]

# extract the atom names for the basis set.
gaussian_atoms = " ".join(set(getAxis(unsorted_geometry,0))) + ' 0'

# generate the input files for nics(1)zz
for i in range(0,len(ssr)):
	# define the name of the files to be generated.
	file_name = str(input)
	file_name = file_name[0:-4]
	file_name = file_name + '_ring' + str(i + 1) + '.com'
	# print(file_name)
	ring_elements = list(ssr[i])
	ring_elements = [i + 1 for i in ring_elements]
	ring_elements = " ".join(map(str, ring_elements))
	# insert the NICS probe and write the file
	nics = insertNICS(geometry,ring_elements)
	orig_stdout = sys.stdout
	sys.stdout = open(file_name,'wt')
	print(gaussian_card)
	print()
	print('This is a comment line.')
	print()
	print(charge_and_multiplicity)
	print2D(nics,False)
	print()
	print(gaussian_atoms)
	print(basis)
	print('****')
	print()
	sys.stdout.close()
	sys.stdout=orig_stdout
	print(ring_elements)
	print2D(nics, False)


delete_file('temp.xyz')
delete_file('temp.mol')