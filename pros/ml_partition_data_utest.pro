pro ml_partition_data_utest

; Test 1
vsini = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
ntraining = 5

sets = ml_partition_data(vsini, ntraining)

; Set2 should be 8, 3, 6
; Set1 should be 9, 2
; Set0 should be 10, 1, 7, 5, 4
s2 = [8,3, 6]
s1 = [9,2]
s0 = [10,1,7,5,4]

printcol, s0, vsini[sets.set0]
print, ''
printcol, s1, vsini[sets.set1]
print, ''
printcol, s2, vsini[sets.set2]
print, ''

stop

; Test 2
data = [[10,9,8,7,6,5,4,3,2,1], $
        [1,2,3,4,5,6,7,8,9,10]]

ntraining = [4, 4, 2]

sets = ml_partition_data(data, ntraining)

s0 = [0,9,6,3]
s1 = [1,8,5,4]
s2 = [2,7]

printcol, s0, sets.set0
print, ''
printcol, s1, sets.set1
print, ''
printcol, s2, sets.set2
print, ''

stop

; Test 3
data = [[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], $
        [1,1,1,5,5,5,2,2,2,21,20, 4,25, 1,19], $
        [1,1,1,5,5,5,2,2,2, 5,20, 4,25, 1,19]]

ntraining = [5, 4, 3, 3]

sets = ml_partition_data(data, ntraining)

s0 = [14,0,9,6,8]
s1 = [13,1,10,7]
s2 = [12,2,4]
s3 = [11,3,5]

printcol, s0, sets.set0
print, ''
printcol, s1, sets.set1
print, ''
printcol, s2, sets.set2
print, ''
printcol, s3, sets.set3

stop

; Test 4
data = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

ntraining = [5,5,5]

sets = ml_partition_data(data, ntraining, /span)

s0 = [14,0]
s1 = [13,1]
s2 = [12,2]

printcol, s0, sets.set0[0:1]
print, ''
printcol, s1, sets.set1[0:1]
print, ''
printcol, s2, sets.set2[0:1]
print, ''
printcol, sets.set0, sets.set1, sets.set2

stop

end
