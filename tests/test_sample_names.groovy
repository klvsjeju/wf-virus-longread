def naming = new GroovyClassLoader().parseClass(new File('lib/SampleNames.groovy'))

assert naming.singleName(null, 'isolate') == 'isolate'
assert naming.singleName('Chosen Name', 'isolate') == 'Chosen_Name'

def folders = [
    ['barcode01', ['one.fastq']],
    ['Patient Two', ['two.fastq']],
    ['extra', ['extra.fastq']]
]
def sheet = [
    [barcode: 'barcode01', alias: 'Patient One'],
    [barcode: 'barcode02', alias: 'Patient Two'],
    [barcode: 'barcode03', alias: 'Missing Patient']
]
def all = naming.selectDirectory(folders, sheet, null)
assert all.samples*.getAt(0) == ['Patient_One', 'Patient_Two']
assert all.missing == [[barcode: 'barcode03', alias: 'Missing_Patient', status: 'missing_input']]
assert all.warnings.any { it.contains('extra') }
assert naming.selectDirectory(folders, sheet, 'barcode02').samples == [['Patient_Two', ['two.fastq']]]
assert naming.selectDirectory(folders, sheet, 'Patient Two').samples == [['Patient_Two', ['two.fastq']]]
assert naming.selectDirectory(folders, null, 'barcode01').samples == [['barcode01', ['one.fastq']]]
assert naming.selectDirectory(folders, null, null).samples*.getAt(0) == ['barcode01', 'Patient_Two', 'extra']
assert naming.selectDirectory(folders, sheet, 'Missing Patient').missing.size() == 1

def rejects = { Closure action, String message ->
    try {
        action()
        assert false: "Expected rejection containing ${message}"
    } catch (IllegalArgumentException error) {
        assert error.message.contains(message)
    }
}
rejects({ naming.selectDirectory(folders + [['Patient One', ['conflict.fastq']]], sheet, null) }, 'Conflicting folders')
rejects({ naming.validateSheet(sheet + [sheet[0]]) }, 'unique')
rejects({ naming.validateSheet([[barcode:'barcode01', alias:'barcode03']]) }, 'must not start')
rejects({ naming.validateSheet([[barcode:'sample01', alias:'Patient']]) }, 'at least two digits')
rejects({ naming.validateSheet([[barcode:'barcode01', alias:'Patient', type:'invalid']]) }, 'invalid sample type')
println 'Sample naming assertions passed'
