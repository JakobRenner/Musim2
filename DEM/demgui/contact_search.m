function [ number_contacts ] = contact_search(filename, radius)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pos = csvread(filename);
counter = 0
pos_size = size(pos,1)
for i = 1:pos_size
    for j = i+1:pos_size
        distance = norm( pos(i,:)-pos(j,:) );
        if distance < radius*2
            counter = counter + 1
        end
    end
number_contacts = counter
end


end