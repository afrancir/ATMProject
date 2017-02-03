function logical = isKey1GreaterKey2(key1,key2)
% Description: it compares keys in lexicographical order as explained in
% the paper. if key1>key2 logical=1, otherwise if key1<key2 logical=0.
%
% INPUT : key_i are 1x2 vectors. 
% OUTPUT: logical 1 or 0.


if key1(1)>key2(1)
    logical = 1;
elseif (key1(1)==key2(1)) && (key1(2)>key2(2))
    logical = 1;
% elseif (key1(1)==key2(1)) && (key1(2)==key2(2)) %error check
    % display('Error at isKey1GreaterKey2: both keys have the same values ')
else
    logical = 0;
end

end 
