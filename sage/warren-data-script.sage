# coding=utf-8

# Some "global" variables we'll need
level = 0
ranks = {}

# Fields from labels
fields = {}

def field_data(s):
	r"""
	Returns full field data from field label.
	"""
	deg, r1, abs_disc, n = [int(c) for c in s.split(".")]
	sig = [r1, (deg-r1)//2]
	return [s, deg, sig, abs_disc]

def field_from_label(lab):
	if lab in fields:
		return fields[lab]
	dummy, deg, sig, abs_disc = field_data(lab)
	d = ZZ(abs_disc)
	if sig[0]==0: d=-d
	x = polygen(QQ)
	t = d%4
	assert t in [0,1]
	pol = x^2 - t*x + (t-d)/4
	K = NumberField(pol, 'a')
	fields[lab] = K
	return K

def extra_info(infile, outfile):
	f = open(infile, "r")
	g = open(outfile, "w")
	for line in f:
		# Parse file
		field, label, isog_letter, isog_number, conductor_label, norm, a1, a2, a3, a4, a6, disc, q_curve = line.split()
		
		# Get field
		k = field_from_label(field)
		if field not in ranks:
			ranks[field] = {}
		
		# Get rank & gens / rank bounds for the curve
		rb = "?"
		gens = "0"
		sha = "?"
		
		# Magma
		magma.eval("QQ<x> := PolynomialRing(RationalField());")
		magma.eval("k<a> := NumberField(%s);" % k.polynomial())
		magma.eval("E := BaseChange(EllipticCurve([%s,%s,%s,%s,%s]), k);" % tuple([int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), [a1, a2, a3, a4, a6])]))
		
		# Analytic Rank
		magma.eval("an := AnalyticRank(E);")
		an = int(round(float(magma.eval("an;"))))
		
		isog_class = "%s-%s%s" % (label, isog_letter, isog_number)
		try:
			isog_rank = ranks[field][isog_class]
		except KeyError:
			isog_rank = None
		if isog_rank == "0":
			r = "0"
			rb = "[0,0]"
			
			# Sha
			magma.eval("sha := ConjecturalSha(E, []);")
			sha = str(int(round(float(magma.eval("sha;")))))
		else:
			magma.eval("rb, gens := DescentInformation(E : RankOnly := true, Silent := true);");
			
			lb = magma.eval("rb[1]")
			ub = magma.eval("rb[2]")
			
			if lb != ub:
				if isog_rank != None:
					r = isog_rank
					rb = "[%s,%s]" % (isog_rank, isog_rank)
				else:
					rb = "[%s,%s]" % (lb, ub)
					r = "?"
			elif lb == "0":
				r = "0"
				gens = "0"
				rb = "[0,0]"
				ranks[field][isog_class] = "0"
				
				# Sha
				magma.eval("sha := ConjecturalSha(E, []);")
				sha = str(int(round(float(magma.eval("sha;")))))
			else:
				r = lb
				ranks[field][isog_class] = r
				rb = "[%s,%s]" % (lb, ub)
				ngens = magma.eval("#gens;");
				if ngens != "0":
					list_of_gens = []
					for j in xrange(1, int(ngens) + 1):
						list_of_gens.append("[%s,%s:%s,%s:%s,%s]" % (magma.eval("gens[%i][1][1];" % j), magma.eval("gens[%i][1][2];" % j), magma.eval("gens[%i][2][1];" % j), magma.eval("gens[%i][2][2];" % j), magma.eval("gens[%i][3][1];" % j), magma.eval("gens[%i][3][2];" % j)))
					gens = ngens + " " + " ".join(list_of_gens)
				
				# Sha
				if ngens == r:
					magma.eval("sha := ConjecturalSha(E, gens);")
					sha = str(int(round(float(magma.eval("sha;")))))
		
		# Print out the curve
		print "%s %s %s %s %s %s %s %s %s" % (field, label, isog_letter, isog_number, r, rb, an, gens, sha)
		g.write("%s %s %s %s %s %s %s %s %s\n" % (field, label, isog_letter, isog_number, r, rb, an, gens, sha))
		g.flush()
	g.close()
	f.close()