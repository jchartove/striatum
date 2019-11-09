%obtain the files with .fig extension
files = dir('*.fig');

%loop through the .fig files
for i=1:length(files)

   %obtain the filename of the .fig file
   filename = files(i).name;

   %open the figure without plotting it
   fig = openfig(filename, 'new', 'invisible');

   %save the figure as a jpg
   saveas(fig, strcat(files(i).name(1:end-4), '.png'));

   %close the figure so that the next could be opened without some java problem
   close;

end