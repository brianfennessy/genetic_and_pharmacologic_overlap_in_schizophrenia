# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES NEEDED:
# 
# 1. 'idList'
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILES CREATED:
# 
# 1. Individual files with information on 11292 drugbank IDs
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# idlist is just a list of all drugbank DB IDs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ ME
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this code maps drugbank IDs to one another, reducing overall number of DB IDs to map to other databases
## EXAMPLE: human secretin
##  -- https://www.drugbank.ca/drugs/DB09532
##  -- on this page says this ID (DB09532) also maps to other IDSs (BTD00039, BIOD00039, DB06369, DB00021)
##  -- page for 1 of these other IDs (DB00021) automatically redirects to page for DB09532
##  -- yet, in drugbank data files (eg, full_database.xml), these IDs are not clearly mapped to each other
##  -- we emailed drugbank, they said they should be mapped in the full_database.xml file; alas, they arent
##  -- this code finds a workaround

cd /sc/arion/projects/psychgen/methods/rx/files/drugbank_id_harmonization

for i in `cat idlist`
do wget https://www.drugbank.ca/drugs/${i}
done

ls DB* | wc -l #11292
wc -l idlist #11292

ls DB* | sort | uniq > tmp1
sort idlist | uniq > tmp2

comm -12 tmp1 tmp2 | wc -l  #11292
