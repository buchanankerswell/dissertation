-- don't do anything unless we target latex
if FORMAT ~= "latex" then
  return {}
end

local latex_figure_start = [[
\begin{figure}
\centering
\includegraphics{%s}
]]
local latex_figure_end = '\\end{figure}'

function latex(str)
  return pandoc.RawInline('latex', str)
end

function make_caption(long_caption, short_caption)
  local caption = short_caption
  table.insert(caption, 1, latex('\\caption['))
  table.insert(caption, latex(']{'))
  for i = 1, #long_caption do
    caption[#caption + 1] = long_caption[i]
  end
  table.insert(caption, latex('}\n'))
  return caption
end

function is_image_with_title(img)
  return img.t == "Image" and img.title
end

function Para(para)
  local img = para.content[1]
  if not (#para.content == 1 and is_image_with_title(img)) then
    return nil
  end
  local title = img.title:gsub('^fig:', '')
  local title_inlines = pandoc.read(title).blocks[1].content
  local figure = make_caption(img.caption, title_inlines)
  local fig_start = latex(latex_figure_start:format(img.src))
  local fig_end = latex(latex_figure_end)
  table.insert(figure, 1, fig_start)
  table.insert(figure, fig_end)
  return pandoc.Plain(figure)
end