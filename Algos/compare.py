import sys

# get args from commandline
filename_1 = sys.argv[1] # references 
filename_2 = sys.argv[2] # reconstructed file
strand_len = int(sys.argv[3])

# Open File in Read Mode
file_1 = open(filename_1, 'r')
file_2 = open(filename_2, 'r')
  
file_1_line = file_1.readline()
file_2_line = file_2.readline()
  
# Use as a Counter
line_no = 0
count_char_diff = 0

while file_1_line != '' or file_2_line != '':
  
    # Removing whitespaces
    file_1_line = file_1_line.rstrip()
    file_2_line = file_2_line.rstrip()
  
    # Compare the lines from both file
    if file_1_line != file_2_line:
        count_char = sum(1 for a, b in zip(file_1_line, file_2_line) if a != b) + abs(len(file_1_line) - len(file_1_line))
        count_char_diff = count_char_diff + count_char
  
    # Read the next line from the file
    file_1_line = file_1.readline()
    file_2_line = file_2.readline()
  
    line_no += 1
  
file_1.close()
file_2.close()

with open(filename_1) as file1:
    with open(filename_2) as file2:
        same = set(file1).intersection(file2)
    
count_same_line = len(same) 

print('Per strand accuracy = ' + str((count_same_line / line_no) * 100) + '%\n')


same_char = (line_no * strand_len) - count_char_diff
char_accuracy = (same_char / ( line_no * strand_len)) * 100
print("Char accuracy = "+ str(char_accuracy))
