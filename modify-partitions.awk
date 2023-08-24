#!/usr/bin/awk -f

# Read the gene mapping file into an array
FNR == NR {
    split($2, uces, ",");
    for (i in uces) {
        gene_map[uces[i]] = $1;
    }
    next;
}

# Process the partition file
/^charset/ {
    charset_name = $2;
    charset_ranges = substr($0, index($0, "=") + 2);
    if (charset_name in gene_map) {
        gene = gene_map[charset_name];
        if (gene in new_charset) {
            new_charset[gene] = new_charset[gene] " " charset_ranges;
        } else {
            new_charset[gene] = "charset " gene " = " charset_ranges;
        }
    } else {
        new_charset[charset_name] = $0;
    }
}

END {
    print "#nexus";
    print "begin sets;";
    for (gene in new_charset) {
        gsub(";", "", new_charset[gene]);
        print new_charset[gene] ";";
    }
    print "done;";
}
