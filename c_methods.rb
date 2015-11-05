require "ffi"

class MismatchInfo < FFI::Struct
  layout :num_windows, :int,
         :xposns, :pointer,
         :mismatches, :pointer
end

module CMethods
  extend FFI::Library
  ffi_lib "#{File.dirname(__FILE__)}/lib/methods.so"

  # attach_function :iupac_match, [:string, :string], :int
  attach_function :perc_mismatch, [:string, :string], :double
  attach_function :windowed_str_mismatch,
                  [:string, :string],
                  MismatchInfo.by_ref

  def self.windowed_str_mismatch_wrapper str1, str2
    minfo = CMethods.windowed_str_mismatch str1, str2

    num = minfo[:num_windows]

    [minfo[:xposns].read_array_of_double(num),
     minfo[:mismatches].read_array_of_double(num)].transpose
  end
end
