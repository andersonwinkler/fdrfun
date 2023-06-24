function y = isoctave()
persistent isoct;
if isempty(isoct)
    isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end
y = isoct;