function [presence] = hubnodes_color_size(hubnodes1,hubnodes2)

% 0 = not hub anywhere
% 1 = common hub
% 2 = hub only in hubnodes1 

node_labels = get_label_nodes("AAL116_labels.txt");
presence=zeros(size(node_labels))';

for i=1:length(node_labels)
    if ismember(node_labels(i),hubnodes1)
        presence(i)=presence(i)+2;
        if ismember(node_labels(i),hubnodes2)
            presence(i)=presence(i)-1;
        end
    end
end
presence=[presence presence];

end
