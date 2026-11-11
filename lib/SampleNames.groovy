/** Sample naming and folder matching adapted from wf-bacterial-genomes v2.0.2 ingress. */
class SampleNames {
    static String normalize(String alias) {
        alias.replace(' ', '_')
    }

    static String singleName(String requested, String derived) {
        normalize(requested ?: derived)
    }

    static void validateSheet(List records) {
        if (!records || !records[0].containsKey('barcode') || !records[0].containsKey('alias')) {
            throw new IllegalArgumentException('CSV needs barcode,alias columns and at least one row')
        }
        def barcodes = records.collect { it.barcode }
        def aliases = records.collect { it.alias }
        if (barcodes.any { !(it ==~ /barcode\d\d+/) }) {
            throw new IllegalArgumentException('CSV barcodes must use barcode followed by at least two digits')
        }
        if (barcodes.collect { it.size() }.toSet().size() > 1) {
            throw new IllegalArgumentException('CSV barcodes must all have the same length')
        }
        if (aliases.any { !it || it.startsWith('barcode') }) {
            throw new IllegalArgumentException('CSV aliases must be nonempty and must not start with barcode')
        }
        if (barcodes.toSet().size() != barcodes.size() || aliases.toSet().size() != aliases.size()) {
            throw new IllegalArgumentException('CSV barcode and alias values must be unique')
        }
        def allowedTypes = ['test_sample', 'positive_control', 'negative_control', 'no_template_control']
        if (records[0].containsKey('type') && records.any { !(it.type in allowedTypes) }) {
            throw new IllegalArgumentException('CSV contains an invalid sample type')
        }
        if (records[0].containsKey('analysis_group') && records.any { !it.analysis_group }) {
            throw new IllegalArgumentException('CSV analysis_group needs a value in every row')
        }
    }

    /** folders contains [original folder baseName, read paths] pairs. */
    static Map selectDirectory(List folders, List records, String requested) {
        def selected = []
        def missing = []
        def warnings = []
        def metadata = [:]
        if (records == null) {
            folders.each { folder ->
                if (!requested || folder[0] == requested) {
                    selected << [normalize(folder[0]), folder[1]]
                    metadata[normalize(folder[0])] = [barcode: folder[0]]
                } else {
                    warnings << "Ignoring ${folder[0]}: folder does not match sample ${requested}"
                }
            }
        } else {
            validateSheet(records)
            def names = folders.collect { it[0] }
            if (names.toSet().size() != names.size()) {
                throw new IllegalArgumentException('Input folders have ambiguous base names')
            }
            def byName = folders.collectEntries { [(it[0]): it[1]] }
            def known = records.collectMany { [it.barcode, it.alias] }.toSet()
            folders.findAll { !(it[0] in known) }.each {
                warnings << "Ignoring ${it[0]}: folder has no sample sheet entry"
            }
            records.each { row ->
                def barcodeReads = byName[row.barcode]
                def aliasReads = byName[row.alias]
                if (barcodeReads && aliasReads) {
                    throw new IllegalArgumentException(
                        "Conflicting folders ${row.barcode} and ${row.alias} for the same sample")
                }
                def reads = barcodeReads ?: aliasReads
                if (!reads) {
                    warnings << "Ignoring ${row.alias}: sample sheet entry has no FASTQ folder"
                }
                if (!requested || requested == row.alias || requested == row.barcode) {
                    if (reads) {
                        selected << [normalize(row.alias), reads]
                        metadata[normalize(row.alias)] = row
                    } else {
                        missing << [barcode: row.barcode, alias: normalize(row.alias), status: 'missing_input']
                    }
                } else if (reads) {
                    warnings << "Ignoring ${row.alias}: sample sheet entry does not match sample ${requested}"
                }
            }
        }
        [samples: selected, missing: missing, warnings: warnings, metadata: metadata]
    }
}
