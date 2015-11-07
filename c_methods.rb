require "ffi"

class MismatchInfo < FFI::Struct
  layout :num_windows, :int,
         :xposns, :pointer,
         :mismatches, :pointer
end

# module LibC
#   extend FFI::Library
#   ffi_lib FFI::Library::LIBC
#   attach_function :malloc, [:size_t], :pointer

#   attach_function :malloc, [:size_t], :pointer
#   attach_function :free, [:pointer], :void
# end

module CMethods
  extend FFI::Library
  ffi_lib "#{File.dirname(__FILE__)}/lib/methods.so"

  # attach_function :iupac_match, [:string, :string], :int
  attach_function :perc_mismatch, [:string, :string], :double
  attach_function :windowed_str_mismatch,
                  [:string, :string],
                  MismatchInfo.by_ref

  attach_function :get_de,
                  [:pointer, :pointer, :int],
                  :double

  def self.get_de_wrapper obs_perc_diffs, exp_perc_diffs
    opd_vals = obs_perc_diffs.map(&:last)
    epd_vals = exp_perc_diffs.map(&:last)

    assert opd_vals.count == epd_vals.count
    opd_pointer = FFI::MemoryPointer.new(:double, opd_vals.count)
    opd_pointer.write_array_of_double opd_vals

    epd_pointer = FFI::MemoryPointer.new(:double, opd_vals.count)
    epd_pointer.write_array_of_double epd_vals

    self.get_de opd_pointer, epd_pointer, opd_vals.count
  end

  def self.windowed_str_mismatch_wrapper str1, str2
    minfo = CMethods.windowed_str_mismatch str1, str2

    num = minfo[:num_windows]

    [minfo[:xposns].read_array_of_double(num),
     minfo[:mismatches].read_array_of_double(num)].transpose
  end
end
