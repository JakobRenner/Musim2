function [ number_contacts ] = contact_search(filename, radius)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pos = csvread(filename);
counter = 0
for i = 1:size(pos,1)-1
    distance = norm( pos(i,:)-pos(i+1,:) ):
    size(distance,1)
    if distance(:) < radius*2
        counter = counter + 1
    end
end
number_contacts = counter
end
