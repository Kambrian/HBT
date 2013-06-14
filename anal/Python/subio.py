class Catalog:
	def __init__(self):
		pass

class GroupCatalog(Catalog):
	def __init__(self):
		Catalog.__init__(self)
		GroupCatalog.Ngroups=0
		
class SubCatalog(Catalog):
