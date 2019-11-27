import csv
import geopandas as gpd

shape_file = '/home/steinar/Documents/edinburgh/data/naturalearth/ne_10m_admin_0_countries.shp'
#continent = 'europe'


world = gpd.read_file(shape_file)
world_info = [iso for iso in world["geometry"]]
world_names = [iso for iso in world["NAME"]]
polygon1 = world_info[0]
polygon2 = world_info[1]
polygon187 = world_info[187]

print(world_names[0], world_names[187])

#print(polygon1.distance(polygon187), polygon2)


world_isos = [iso for iso in world["ISO_A3"]]

ISO_to_population = {}
ISO_to_GDP = {}
ISO_to_POLYGON = {}

for ISO, POP_EST, GDP, POLYGON in zip(world["ISO_A3"], world["POP_EST"], world["GDP_MD_EST"], world["geometry"]):
    ISO_to_population[ISO] = int(POP_EST)
    ISO_to_GDP[ISO] = int(GDP)
    ISO_to_POLYGON[ISO] = POLYGON

writer = csv.writer(open('pairwise_distances.csv', mode='w'))

for iso1 in world_isos:
    for iso2 in world_isos:
        if iso1 != '-99' and iso2 != '-99':
            print(iso1, iso2)
            reporter_poly = ISO_to_POLYGON[iso1]
            partner_poly = ISO_to_POLYGON[iso2]
            distance = reporter_poly.distance(partner_poly)
            print([iso1, iso2, distance])
            writer.writerow([iso1, iso2, distance])
            



            
