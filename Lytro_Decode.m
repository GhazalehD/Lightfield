function [  ] = Lytro_Decode( Path )
DecodeOptions.OptionalTasks = 'ColourCorrect';
LFUtilDecodeLytroFolder(Path, [], DecodeOptions);
end

