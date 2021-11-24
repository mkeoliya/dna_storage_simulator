import noisy_strand_gen as gen
#suppose we want to transform string2 into string1
#first return value is the edit distance, second return value is the list of operations
def edit_distance(s1, s2, m, n, matrix, matrix2):
 
    if matrix[m][n] != -1:
        return matrix[m][n], matrix2[m][n]
    if m == 0:
        matrix[m][n] = n
        ops = ['-' for x in range(n)]
        matrix2[m][n] = ops
        return n, ops
    if n == 0:
        matrix[m][n] = m
        ops = ['+' for x in range(m)]
        matrix2[m][n] = ops 
        return m, ops
    if s1[m-1] == s2[n-1]:
        dist, ops = edit_distance(s1[:m-1], s2[:n-1], m-1, n-1, matrix, matrix2)
        ops = ops.copy()
        ops.append('=')
        matrix[m][n] = dist
        matrix2[m][n] = ops
        return dist, ops
    else:
        dist1, ops1 = edit_distance(s1[:m-1], s2[:n], m-1, n, matrix, matrix2)
        dist2, ops2 = edit_distance(s1[:m], s2[:n-1], m, n-1, matrix, matrix2)
        dist3, ops3 = edit_distance(s1[:m-1], s2[:n-1], m-1, n-1, matrix, matrix2)
        dist = min(dist1, dist2, dist3)
        if dist == dist1:
            ops1 = ops1.copy()
            ops1.append('+')
            matrix[m][n] = dist1 + 1
            matrix2[m][n] = ops1
      
            return dist1 + 1, ops1
        elif dist == dist2:
            ops2 = ops2.copy()
            ops2.append('-')
            matrix[m][n] = dist2 + 1
            matrix2[m][n] = ops2
     
            return dist2 + 1, ops2
        else:
            ops3 = ops3.copy()
            ops3.append('s')
            matrix[m][n] = dist3 + 1
            matrix2[m][n] = ops3
      
            return dist3 + 1, ops3

def get_reads(reads_file):
    return
def get_references(references_file):
    return
def match(reads, references):
    read_matches = []
    ref_matches = []
    for i in range(len(references)):
        ref_matches.append([])
    ops_to_transform = []
    for i in range(len(reads)):
        read = reads[i]
        best_ref_index = -1
        lowest_dist = 1000000
        best_ops = []
        for j in range(len(references)):
            ref = references[j]
            m = len(read)
            n = len(ref)
            matrix = [[-1 for x in range(n+1)] for y in range(m+1)]
            matrix2 = [[-1 for x in range(n+1)] for y in range(m+1)]
            dist, ops = edit_distance(read, ref, m, n, matrix, matrix2)
            if dist < lowest_dist:
                lowest_dist = dist
                best_ref_index = j
                best_ops = ops
        read_matches.append(best_ref_index)
        ref_matches[best_ref_index].append((i, lowest_dist))
        ops_to_transform.append(best_ops)
    return read_matches, ref_matches, ops_to_transform
    
# returns 2 items, (edit distance, list of operations to transform s2 into s1)
def ops_list(s1, s2):
    m = len(s1)
    n = len(s2)
    matrix = [[-1 for x in range(n+1)] for y in range(m+1)]
    matrix2 = [[-1 for x in range(n+1)] for y in range(m+1)]
    return edit_distance(s1, s2, m, n, matrix, matrix2)

    
    
    
    
