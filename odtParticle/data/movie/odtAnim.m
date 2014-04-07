% D. O. Lignell, 2/19/09
% Create a movie of ODT data

clc;
clear;

myfile = 'fileList';             % list of file names to read

f1 = fopen(myfile);          

i = 1;
while(1);
    file = fscanf(f1, '%s', 1);  % get the name of each file
    if(feof(f1)) break; end
        flist{i,1} = file;
        i = i+1;
end
fclose(f1);                      % close the file name list

nfiles = length(flist);          % number of files to read

box('on');

for ifi=1:nfiles/8                 % read each file and get its data
    ifi
    f1 = fopen(flist{ifi});
    ln = fgetl(f1);              % dummy header lines
    ln = fgetl(f1);              % dummy header lines 
    titles = fgetl(f1);          % names of columns
    i=1;
    while(~feof(f1))             % get data in each row, fill matrix A
        ln = fgetl(f1);
        A(i,:) = [sscanf(ln,'%f')]';
        i = i+1;
    end
    fclose(f1);                  % close file
    %plot(A(:,6), A(:,1),'.-');    % plot profile
    %axis([0 22 0 2]);            % set axis
    plot(A(:,1), A(:,6),'.-','LineWidth', 1.5);    % plot profile
    

    axis([0 2 0 22]);            % set axis
    M(ifi) = getframe;           % store the picture
    clear A;                     % clear the data array
end

%movie(M);                        % replay the movie
map = colormap;
mpgwrite(M,map,'myanim.mpg');    % create the external movie






