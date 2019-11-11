'''
Extracts the final image from an <.xyz> trajectory file.

'''

file_name = str(input('Name of trajectory.xyz file: '))

image = []
with open(file_name, 'r') as f:
    for n, line in enumerate(f):
        if (' i =') in line:
            image_line_number = n

with open(file_name, 'r') as f:
    for n, line in enumerate(f):
        if n >= image_line_number:
            image.append(line)

with open('FINXYZ.xyz', 'a+') as f:
    f.write('{}\n'.format('Final Image'))
    for i in range(0, len(image)):
        f.write(str(image[i]))
        
