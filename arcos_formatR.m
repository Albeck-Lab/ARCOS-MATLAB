function arcos_formatR(filename,XCoord,YCoord,bin)
    array = zeros(size(bin,1)*size(bin,2),6);
    sz = size(bin,1);
    for col = 1:size(bin,2)
        for row = 1:size(bin,1)
            if col >1
                ind = row+sz*(col-1); 
            else
                ind = row;
            end
            array(ind,1) = col;
            array(ind,2) = XCoord(row,col);
            array(ind,3) = YCoord(row,col);
            array(ind,4) = bin(row,col);
            array(ind,5) = row;
        end
    end
    table = array2table(array);
    table.Properties.VariableNames(1:5) = ["t","x","y","m","id"];
    writetable(table,append(filename,'.csv'));
end