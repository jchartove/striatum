function datatable = spectramean(structure,filenums)
    datatable = zeros(length(filenums),151);
    iterator = 0;
    for num = filenums
        iterator = iterator + 1;
        datatable(iterator,:) = structure(num).values;
    end
end