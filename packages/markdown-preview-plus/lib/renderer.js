"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const path = require("path");
const pandocHelper = require("./pandoc-helper");
const markdownIt = require("./markdown-it-helper");
const extension_helper_1 = require("./extension-helper");
const atom_1 = require("atom");
const util_1 = require("./util");
const util_common_1 = require("./util-common");
const { resourcePath } = atom.getLoadSettings();
const packagePath = path.dirname(__dirname);
async function render(options) {
    const text = options.text.replace(/^\s*<!doctype(\s+.*)?>\s*/i, '');
    let html;
    let error;
    if (util_1.atomConfig().renderer === 'pandoc') {
        try {
            html = await pandocHelper.renderPandoc(text, options.filePath, options.renderLaTeX);
        }
        catch (err) {
            const e = err;
            if (e.html === undefined)
                throw e;
            error = e.message;
            html = e.html;
        }
    }
    else {
        html = markdownIt.render(text, options.renderLaTeX);
    }
    const parser = new DOMParser();
    const doc = parser.parseFromString(html, 'text/html');
    sanitize(doc);
    if (options.mode === 'normal') {
        if (options.imageWatcher)
            options.imageWatcher.clear();
        resolveImagePaths(doc, options.filePath, false, undefined, options.imageWatcher);
    }
    else {
        switch (options.mode) {
            case 'save':
                handleImages({
                    doc,
                    filePath: options.filePath,
                    savePath: options.savePath,
                    behaviour: util_1.atomConfig().saveConfig.mediaOnSaveAsHTMLBehaviour,
                });
                break;
            case 'copy':
                handleImages({
                    doc,
                    filePath: options.filePath,
                    behaviour: util_1.atomConfig().saveConfig.mediaOnCopyAsHTMLBehaviour,
                });
                break;
            default:
                throw invalidMode(options);
        }
    }
    let defaultCodeLanguage = 'text';
    if ((options.grammar && options.grammar.scopeName) === 'source.litcoffee') {
        defaultCodeLanguage = 'coffee';
    }
    if (!(util_1.atomConfig().renderer === 'pandoc' &&
        util_1.atomConfig().pandocConfig.useNativePandocCodeStyles)) {
        await highlightCodeBlocks(doc, defaultCodeLanguage);
    }
    if (error) {
        const errd = doc.createElement('div');
        const msgel = doc.createElement('code');
        msgel.innerText = error;
        errd.innerHTML = `<h1>Pandoc Error:</h1>${msgel.outerHTML}<hr>`;
        doc.body.insertBefore(errd, doc.body.firstElementChild);
    }
    return doc;
}
exports.render = render;
function invalidMode(mode) {
    return new Error(`Invalid render mode ${JSON.stringify(mode)}`);
}
function sanitize(doc) {
    doc.querySelectorAll("script:not([type^='math/tex'])").forEach((elem) => {
        elem.remove();
    });
    const attributesToRemove = [
        'onabort',
        'onblur',
        'onchange',
        'onclick',
        'ondbclick',
        'onerror',
        'onfocus',
        'onkeydown',
        'onkeypress',
        'onkeyup',
        'onload',
        'onmousedown',
        'onmousemove',
        'onmouseover',
        'onmouseout',
        'onmouseup',
        'onreset',
        'onresize',
        'onscroll',
        'onselect',
        'onsubmit',
        'onunload',
    ];
    doc.querySelectorAll('*').forEach((elem) => attributesToRemove.map((attribute) => {
        elem.removeAttribute(attribute);
    }));
}
function handleImages(opts) {
    const relativize = opts.behaviour === 'relativized';
    switch (opts.behaviour) {
        case 'relativized':
        case 'absolutized':
            resolveImagePaths(opts.doc, opts.filePath, relativize, opts.savePath);
            break;
        case 'untouched':
    }
}
function resolveImagePaths(doc, filePath, relativize, savePath, imageWatcher) {
    const [rootDirectory] = atom.project.relativizePath(filePath || '');
    const media = util_common_1.getMedia(doc);
    Array.from(media).map(function (img) {
        let attrName;
        if (img.tagName === 'LINK')
            attrName = 'href';
        else
            attrName = 'src';
        let src = img.getAttribute(attrName);
        if (src) {
            if (util_1.atomConfig().renderer !== 'pandoc') {
                src = decodeURI(src);
            }
            if (src.match(/^(https?|atom|data):/)) {
                return;
            }
            if (process.resourcesPath && src.startsWith(process.resourcesPath)) {
                return;
            }
            if (src.startsWith(resourcePath)) {
                return;
            }
            if (src.startsWith(packagePath)) {
                return;
            }
            if (src[0] === '/') {
                if (!util_1.isFileSync(src)) {
                    try {
                        if (rootDirectory !== null) {
                            src = path.join(rootDirectory, src.substring(1));
                        }
                    }
                    catch (e) {
                    }
                }
            }
            else if (filePath) {
                src = path.resolve(path.dirname(filePath), src);
            }
            if (relativize && (filePath !== undefined || savePath !== undefined)) {
                const fp = savePath !== undefined ? savePath : filePath;
                src = path.relative(path.dirname(fp), src);
            }
            if (imageWatcher) {
                const v = imageWatcher.watch(src);
                if (v !== undefined)
                    src = `${src}?v=${v}`;
            }
            img[attrName] = src;
        }
    });
}
async function highlightCodeBlocks(domFragment, defaultLanguage) {
    const fontFamily = atom.config.get('editor.fontFamily');
    if (fontFamily) {
        for (const codeElement of Array.from(domFragment.querySelectorAll('code'))) {
            codeElement.style.fontFamily = fontFamily;
        }
    }
    await Promise.all(Array.from(domFragment.querySelectorAll('pre')).map(async (preElement) => {
        const codeBlock = preElement.firstElementChild !== null
            ? preElement.firstElementChild
            : preElement;
        const cbClass = codeBlock.className || preElement.className;
        const fenceName = cbClass
            ? cbClass.replace(/^(lang-|sourceCode )/, '')
            : defaultLanguage;
        const ctw = util_1.atomConfig().codeTabWidth;
        const ed = new atom_1.TextEditor({
            readonly: true,
            keyboardInputEnabled: false,
            showInvisibles: false,
            tabLength: ctw === 0 ? atom.config.get('editor.tabLength') : ctw,
        });
        const el = atom.views.getView(ed);
        try {
            el.setUpdatedSynchronously(true);
            el.style.pointerEvents = 'none';
            el.style.position = 'absolute';
            el.style.top = '100vh';
            el.style.width = '100vw';
            atom.grammars.assignLanguageMode(ed.getBuffer(), extension_helper_1.scopeForFenceName(fenceName));
            ed.setText(codeBlock.textContent.replace(/\r?\n$/, ''));
            atom.views.getView(atom.workspace).appendChild(el);
            await editorTokenized(ed);
            const html = Array.from(el.querySelectorAll('.line:not(.dummy)'));
            preElement.classList.add('editor-colors');
            preElement.innerHTML = html.map((x) => x.innerHTML).join('\n');
            if (fenceName)
                preElement.classList.add(`lang-${fenceName}`);
        }
        finally {
            el.remove();
        }
    }));
    return domFragment;
}
async function editorTokenized(editor) {
    return new Promise((resolve) => {
        const languageMode = editor.getBuffer().getLanguageMode();
        const nextUpdatePromise = editor.component.getNextUpdatePromise();
        if (languageMode.fullyTokenized || languageMode.tree) {
            resolve(nextUpdatePromise);
        }
        else {
            const disp = editor.onDidTokenize(() => {
                disp.dispose();
                resolve(nextUpdatePromise);
            });
        }
    });
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoicmVuZGVyZXIuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyIuLi9zcmMvcmVuZGVyZXIudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7QUFBQSw2QkFBNkI7QUFDN0IsZ0RBQWdEO0FBQ2hELG1EQUFtRDtBQUNuRCx5REFBc0Q7QUFDdEQsK0JBQTBDO0FBQzFDLGlDQUErQztBQUMvQywrQ0FBd0M7QUFHeEMsTUFBTSxFQUFFLFlBQVksRUFBRSxHQUFHLElBQUksQ0FBQyxlQUFlLEVBQUUsQ0FBQTtBQUMvQyxNQUFNLFdBQVcsR0FBRyxJQUFJLENBQUMsT0FBTyxDQUFDLFNBQVMsQ0FBQyxDQUFBO0FBaUJwQyxLQUFLLFVBQVUsTUFBTSxDQUFDLE9BQXNCO0lBR2pELE1BQU0sSUFBSSxHQUFHLE9BQU8sQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLDRCQUE0QixFQUFFLEVBQUUsQ0FBQyxDQUFBO0lBRW5FLElBQUksSUFBSSxDQUFBO0lBQ1IsSUFBSSxLQUFLLENBQUE7SUFDVCxJQUFJLGlCQUFVLEVBQUUsQ0FBQyxRQUFRLEtBQUssUUFBUSxFQUFFO1FBQ3RDLElBQUk7WUFDRixJQUFJLEdBQUcsTUFBTSxZQUFZLENBQUMsWUFBWSxDQUNwQyxJQUFJLEVBQ0osT0FBTyxDQUFDLFFBQVEsRUFDaEIsT0FBTyxDQUFDLFdBQVcsQ0FDcEIsQ0FBQTtTQUNGO1FBQUMsT0FBTyxHQUFHLEVBQUU7WUFDWixNQUFNLENBQUMsR0FBRyxHQUFnQyxDQUFBO1lBQzFDLElBQUksQ0FBQyxDQUFDLElBQUksS0FBSyxTQUFTO2dCQUFFLE1BQU0sQ0FBQyxDQUFBO1lBQ2pDLEtBQUssR0FBRyxDQUFDLENBQUMsT0FBaUIsQ0FBQTtZQUMzQixJQUFJLEdBQUcsQ0FBQyxDQUFDLElBQWMsQ0FBQTtTQUN4QjtLQUNGO1NBQU07UUFDTCxJQUFJLEdBQUcsVUFBVSxDQUFDLE1BQU0sQ0FBQyxJQUFJLEVBQUUsT0FBTyxDQUFDLFdBQVcsQ0FBQyxDQUFBO0tBQ3BEO0lBQ0Qs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QixFQUFFLENBQUE7UUFDakUsSUFBSSxZQUFZLENBQUMsY0FBYyxJQUFJLFlBQVksQ0FBQyxJQUFJLEVBQUU7WUFDcEQsT0FBTyxDQUFDLGlCQUFpQixDQUFDLENBQUE7U0FDM0I7YUFBTTtZQUNMLE1BQU0sSUFBSSxHQUFHLE1BQU0sQ0FBQyxhQUFhLENBQUMsR0FBRyxFQUFFO2dCQUNyQyxJQUFJLENBQUMsT0FBTyxFQUFFLENBQUE7Z0JBQ2QsT0FBTyxDQUFDLGlCQUFpQixDQUFDLENBQUE7WUFDNUIsQ0FBQyxDQUFDLENBQUE7U0FDSDtJQUNILENBQUMsQ0FBQyxDQUFBO0FBQ0osQ0FBQyIsInNvdXJjZXNDb250ZW50IjpbImltcG9ydCBwYXRoID0gcmVxdWlyZSgncGF0aCcpXG5pbXBvcnQgcGFuZG9jSGVscGVyID0gcmVxdWlyZSgnLi9wYW5kb2MtaGVscGVyJylcbmltcG9ydCBtYXJrZG93bkl0ID0gcmVxdWlyZSgnLi9tYXJrZG93bi1pdC1oZWxwZXInKSAvLyBEZWZlciB1bnRpbCB1c2VkXG5pbXBvcnQgeyBzY29wZUZvckZlbmNlTmFtZSB9IGZyb20gJy4vZXh0ZW5zaW9uLWhlbHBlcidcbmltcG9ydCB7IEdyYW1tYXIsIFRleHRFZGl0b3IgfSBmcm9tICdhdG9tJ1xuaW1wb3J0IHsgaXNGaWxlU3luYywgYXRvbUNvbmZpZyB9IGZyb20gJy4vdXRpbCdcbmltcG9ydCB7IGdldE1lZGlhIH0gZnJvbSAnLi91dGlsLWNvbW1vbidcbmltcG9ydCB7IEltYWdlV2F0Y2hlciB9IGZyb20gJy4vaW1hZ2Utd2F0Y2gtaGVscGVyJ1xuXG5jb25zdCB7IHJlc291cmNlUGF0aCB9ID0gYXRvbS5nZXRMb2FkU2V0dGluZ3MoKVxuY29uc3QgcGFja2FnZVBhdGggPSBwYXRoLmRpcm5hbWUoX19kaXJuYW1lKVxuXG5leHBvcnQgdHlwZSBSZW5kZXJNb2RlID0gJ25vcm1hbCcgfCAnY29weScgfCAnc2F2ZSdcblxuZXhwb3J0IGludGVyZmFjZSBDb21tb25SZW5kZXJPcHRpb25zPFQgZXh0ZW5kcyBSZW5kZXJNb2RlPiB7XG4gIHRleHQ6IHN0cmluZ1xuICBmaWxlUGF0aDogc3RyaW5nIHwgdW5kZWZpbmVkXG4gIGdyYW1tYXI/OiBHcmFtbWFyXG4gIHJlbmRlckxhVGVYOiBib29sZWFuXG4gIG1vZGU6IFRcbn1cblxuZXhwb3J0IHR5cGUgUmVuZGVyT3B0aW9ucyA9XG4gIHwgKENvbW1vblJlbmRlck9wdGlvbnM8J25vcm1hbCc+ICYgeyBpbWFnZVdhdGNoZXI/OiBJbWFnZVdhdGNoZXIgfSlcbiAgfCAoQ29tbW9uUmVuZGVyT3B0aW9uczwnc2F2ZSc+ICYgeyBzYXZlUGF0aDogc3RyaW5nIH0pXG4gIHwgKENvbW1vblJlbmRlck9wdGlvbnM8J2NvcHknPilcblxuZXhwb3J0IGFzeW5jIGZ1bmN0aW9uIHJlbmRlcihvcHRpb25zOiBSZW5kZXJPcHRpb25zKTogUHJvbWlzZTxIVE1MRG9jdW1lbnQ+IHtcbiAgLy8gUmVtb3ZlIHRoZSA8IWRvY3R5cGU+IHNpbmNlIG90aGVyd2lzZSBtYXJrZWQgd2lsbCBlc2NhcGUgaXRcbiAgLy8gaHR0cHM6Ly9naXRodWIuY29tL2NoamovbWFya2VkL2lzc3Vlcy8zNTRcbiAgY29uc3QgdGV4dCA9IG9wdGlvbnMudGV4dC5yZXBsYWNlKC9eXFxzKjwhZG9jdHlwZShcXHMrLiopPz5cXHMqL2ksICcnKVxuXG4gIGxldCBodG1sXG4gIGxldCBlcnJvclxuICBpZiAoYXRvbUNvbmZpZygpLnJlbmRlcmVyID09PSAncGFuZG9jJykge1xuICAgIHRyeSB7XG4gICAgICBodG1sID0gYXdhaXQgcGFuZG9jSGVscGVyLnJlbmRlclBhbmRvYyhcbiAgICAgICAgdGV4dCxcbiAgICAgICAgb3B0aW9ucy5maWxlUGF0aCxcbiAgICAgICAgb3B0aW9ucy5yZW5kZXJMYVRlWCxcbiAgICAgIClcbiAgICB9IGNhdGNoIChlcnIpIHtcbiAgICAgIGNvbnN0IGUgPSBlcnIgYXMgRXJyb3IgJiB7IGh0bWw/OiBzdHJpbmcgfVxuICAgICAgaWYgKGUuaHRtbCA9PT0gdW5kZWZpbmVkKSB0aHJvdyBlXG4gICAgICBlcnJvciA9IGUubWVzc2FnZSBhcyBzdHJpbmdcbiAgICAgIGh0bWwgPSBlLmh0bWwgYXMgc3RyaW5nXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGh0bWwgPSBtYXJrZG93bkl0LnJlbmRlcih0ZXh0LCBvcHRpb25zLnJlbmRlckxhVGVYKVxuICB9XG4gIGNvbnN0IHBhcnNlciA9IG5ldyBET01QYXJzZXIoKVxuICBjb25zdCBkb2MgPSBwYXJzZXIucGFyc2VGcm9tU3RyaW5nKGh0bWwsICd0ZXh0L2h0bWwnKVxuICBzYW5pdGl6ZShkb2MpXG4gIGlmIChvcHRpb25zLm1vZGUgPT09ICdub3JtYWwnKSB7XG4gICAgaWYgKG9wdGlvbnMuaW1hZ2VXYXRjaGVyKSBvcHRpb25zLmltYWdlV2F0Y2hlci5jbGVhcigpXG4gICAgcmVzb2x2ZUltYWdlUGF0aHMoXG4gICAgICBkb2MsXG4gICAgICBvcHRpb25zLmZpbGVQYXRoLFxuICAgICAgZmFsc2UsXG4gICAgICB1bmRlZmluZWQsXG4gICAgICBvcHRpb25zLmltYWdlV2F0Y2hlcixcbiAgICApXG4gIH0gZWxzZSB7XG4gICAgc3dpdGNoIChvcHRpb25zLm1vZGUpIHtcbiAgICAgIGNhc2UgJ3NhdmUnOlxuICAgICAgICBoYW5kbGVJbWFnZXMoe1xuICAgICAgICAgIGRvYyxcbiAgICAgICAgICBmaWxlUGF0aDogb3B0aW9ucy5maWxlUGF0aCxcbiAgICAgICAgICBzYXZlUGF0aDogb3B0aW9ucy5zYXZlUGF0aCxcbiAgICAgICAgICBiZWhhdmlvdXI6IGF0b21Db25maWcoKS5zYXZlQ29uZmlnLm1lZGlhT25TYXZlQXNIVE1MQmVoYXZpb3VyLFxuICAgICAgICB9KVxuICAgICAgICBicmVha1xuICAgICAgY2FzZSAnY29weSc6XG4gICAgICAgIGhhbmRsZUltYWdlcyh7XG4gICAgICAgICAgZG9jLFxuICAgICAgICAgIGZpbGVQYXRoOiBvcHRpb25zLmZpbGVQYXRoLFxuICAgICAgICAgIGJlaGF2aW91cjogYXRvbUNvbmZpZygpLnNhdmVDb25maWcubWVkaWFPbkNvcHlBc0hUTUxCZWhhdmlvdXIsXG4gICAgICAgIH0pXG4gICAgICAgIGJyZWFrXG4gICAgICBkZWZhdWx0OlxuICAgICAgICB0aHJvdyBpbnZhbGlkTW9kZShvcHRpb25zKVxuICAgIH1cbiAgfVxuICBsZXQgZGVmYXVsdENvZGVMYW5ndWFnZTogc3RyaW5nID0gJ3RleHQnXG4gIC8vIERlZmF1bHQgY29kZSBibG9ja3MgdG8gYmUgY29mZmVlIGluIExpdGVyYXRlIENvZmZlZVNjcmlwdCBmaWxlc1xuICBpZiAoKG9wdGlvbnMuZ3JhbW1hciAmJiBvcHRpb25zLmdyYW1tYXIuc2NvcGVOYW1lKSA9PT0gJ3NvdXJjZS5saXRjb2ZmZWUnKSB7XG4gICAgZGVmYXVsdENvZGVMYW5ndWFnZSA9ICdjb2ZmZWUnXG4gIH1cbiAgaWYgKFxuICAgICEoXG4gICAgICBhdG9tQ29uZmlnKCkucmVuZGVyZXIgPT09ICdwYW5kb2MnICYmXG4gICAgICBhdG9tQ29uZmlnKCkucGFuZG9jQ29uZmlnLnVzZU5hdGl2ZVBhbmRvY0NvZGVTdHlsZXNcbiAgICApXG4gICkge1xuICAgIGF3YWl0IGhpZ2hsaWdodENvZGVCbG9ja3MoZG9jLCBkZWZhdWx0Q29kZUxhbmd1YWdlKVxuICB9XG4gIGlmIChlcnJvcikge1xuICAgIGNvbnN0IGVycmQgPSBkb2MuY3JlYXRlRWxlbWVudCgnZGl2JylcbiAgICBjb25zdCBtc2dlbCA9IGRvYy5jcmVhdGVFbGVtZW50KCdjb2RlJylcbiAgICBtc2dlbC5pbm5lclRleHQgPSBlcnJvclxuICAgIGVycmQuaW5uZXJIVE1MID0gYDxoMT5QYW5kb2MgRXJyb3I6PC9oMT4ke21zZ2VsLm91dGVySFRNTH08aHI+YFxuICAgIGRvYy5ib2R5Lmluc2VydEJlZm9yZShlcnJkLCBkb2MuYm9keS5maXJzdEVsZW1lbnRDaGlsZClcbiAgfVxuICByZXR1cm4gZG9jXG59XG5cbmZ1bmN0aW9uIGludmFsaWRNb2RlKG1vZGU6IG5ldmVyKSB7XG4gIHJldHVybiBuZXcgRXJyb3IoYEludmFsaWQgcmVuZGVyIG1vZGUgJHtKU09OLnN0cmluZ2lmeShtb2RlKX1gKVxufVxuXG5mdW5jdGlvbiBzYW5pdGl6ZShkb2M6IEhUTUxEb2N1bWVudCkge1xuICAvLyBEbyBub3QgcmVtb3ZlIE1hdGhKYXggc2NyaXB0IGRlbGltaXRlZCBibG9ja3NcbiAgZG9jLnF1ZXJ5U2VsZWN0b3JBbGwoXCJzY3JpcHQ6bm90KFt0eXBlXj0nbWF0aC90ZXgnXSlcIikuZm9yRWFjaCgoZWxlbSkgPT4ge1xuICAgIGVsZW0ucmVtb3ZlKClcbiAgfSlcbiAgY29uc3QgYXR0cmlidXRlc1RvUmVtb3ZlID0gW1xuICAgICdvbmFib3J0JyxcbiAgICAnb25ibHVyJyxcbiAgICAnb25jaGFuZ2UnLFxuICAgICdvbmNsaWNrJyxcbiAgICAnb25kYmNsaWNrJyxcbiAgICAnb25lcnJvcicsXG4gICAgJ29uZm9jdXMnLFxuICAgICdvbmtleWRvd24nLFxuICAgICdvbmtleXByZXNzJyxcbiAgICAnb25rZXl1cCcsXG4gICAgJ29ubG9hZCcsXG4gICAgJ29ubW91c2Vkb3duJyxcbiAgICAnb25tb3VzZW1vdmUnLFxuICAgICdvbm1vdXNlb3ZlcicsXG4gICAgJ29ubW91c2VvdXQnLFxuICAgICdvbm1vdXNldXAnLFxuICAgICdvbnJlc2V0JyxcbiAgICAnb25yZXNpemUnLFxuICAgICdvbnNjcm9sbCcsXG4gICAgJ29uc2VsZWN0JyxcbiAgICAnb25zdWJtaXQnLFxuICAgICdvbnVubG9hZCcsXG4gIF1cbiAgZG9jLnF1ZXJ5U2VsZWN0b3JBbGwoJyonKS5mb3JFYWNoKChlbGVtKSA9PlxuICAgIGF0dHJpYnV0ZXNUb1JlbW92ZS5tYXAoKGF0dHJpYnV0ZSkgPT4ge1xuICAgICAgZWxlbS5yZW1vdmVBdHRyaWJ1dGUoYXR0cmlidXRlKVxuICAgIH0pLFxuICApXG59XG5cbmZ1bmN0aW9uIGhhbmRsZUltYWdlcyhvcHRzOiB7XG4gIGJlaGF2aW91cjogJ3JlbGF0aXZpemVkJyB8ICdhYnNvbHV0aXplZCcgfCAndW50b3VjaGVkJ1xuICBkb2M6IEhUTUxEb2N1bWVudFxuICBmaWxlUGF0aD86IHN0cmluZ1xuICBzYXZlUGF0aD86IHN0cmluZ1xufSkge1xuICBjb25zdCByZWxhdGl2aXplID0gb3B0cy5iZWhhdmlvdXIgPT09ICdyZWxhdGl2aXplZCdcbiAgc3dpdGNoIChvcHRzLmJlaGF2aW91cikge1xuICAgIGNhc2UgJ3JlbGF0aXZpemVkJzpcbiAgICBjYXNlICdhYnNvbHV0aXplZCc6XG4gICAgICByZXNvbHZlSW1hZ2VQYXRocyhvcHRzLmRvYywgb3B0cy5maWxlUGF0aCwgcmVsYXRpdml6ZSwgb3B0cy5zYXZlUGF0aClcbiAgICAgIGJyZWFrXG4gICAgY2FzZSAndW50b3VjaGVkJzpcbiAgICAvKiBub29wICovXG4gIH1cbn1cblxuZnVuY3Rpb24gcmVzb2x2ZUltYWdlUGF0aHMoXG4gIGRvYzogSFRNTERvY3VtZW50LFxuICBmaWxlUGF0aDogc3RyaW5nIHwgdW5kZWZpbmVkLFxuICByZWxhdGl2aXplOiBib29sZWFuLFxuICBzYXZlUGF0aD86IHN0cmluZyxcbiAgaW1hZ2VXYXRjaGVyPzogSW1hZ2VXYXRjaGVyLFxuKSB7XG4gIGNvbnN0IFtyb290RGlyZWN0b3J5XSA9IGF0b20ucHJvamVjdC5yZWxhdGl2aXplUGF0aChmaWxlUGF0aCB8fCAnJylcbiAgY29uc3QgbWVkaWEgPSBnZXRNZWRpYShkb2MpXG4gIEFycmF5LmZyb20obWVkaWEpLm1hcChmdW5jdGlvbihpbWcpIHtcbiAgICBsZXQgYXR0ck5hbWU6ICdocmVmJyB8ICdzcmMnXG4gICAgaWYgKGltZy50YWdOYW1lID09PSAnTElOSycpIGF0dHJOYW1lID0gJ2hyZWYnXG4gICAgZWxzZSBhdHRyTmFtZSA9ICdzcmMnXG4gICAgbGV0IHNyYyA9IGltZy5nZXRBdHRyaWJ1dGUoYXR0ck5hbWUpXG4gICAgaWYgKHNyYykge1xuICAgICAgaWYgKGF0b21Db25maWcoKS5yZW5kZXJlciAhPT0gJ3BhbmRvYycpIHtcbiAgICAgICAgc3JjID0gZGVjb2RlVVJJKHNyYylcbiAgICAgIH1cblxuICAgICAgaWYgKHNyYy5tYXRjaCgvXihodHRwcz98YXRvbXxkYXRhKTovKSkge1xuICAgICAgICByZXR1cm5cbiAgICAgIH1cbiAgICAgIGlmIChwcm9jZXNzLnJlc291cmNlc1BhdGggJiYgc3JjLnN0YXJ0c1dpdGgocHJvY2Vzcy5yZXNvdXJjZXNQYXRoKSkge1xuICAgICAgICByZXR1cm5cbiAgICAgIH1cbiAgICAgIGlmIChzcmMuc3RhcnRzV2l0aChyZXNvdXJjZVBhdGgpKSB7XG4gICAgICAgIHJldHVyblxuICAgICAgfVxuICAgICAgaWYgKHNyYy5zdGFydHNXaXRoKHBhY2thZ2VQYXRoKSkge1xuICAgICAgICByZXR1cm5cbiAgICAgIH1cblxuICAgICAgaWYgKHNyY1swXSA9PT0gJy8nKSB7XG4gICAgICAgIGlmICghaXNGaWxlU3luYyhzcmMpKSB7XG4gICAgICAgICAgdHJ5IHtcbiAgICAgICAgICAgIGlmIChyb290RGlyZWN0b3J5ICE9PSBudWxsKSB7XG4gICAgICAgICAgICAgIHNyYyA9IHBhdGguam9pbihyb290RGlyZWN0b3J5LCBzcmMuc3Vic3RyaW5nKDEpKVxuICAgICAgICAgICAgfVxuICAgICAgICAgIH0gY2F0Y2ggKGUpIHtcbiAgICAgICAgICAgIC8vIG5vb3BcbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgIH0gZWxzZSBpZiAoZmlsZVBhdGgpIHtcbiAgICAgICAgc3JjID0gcGF0aC5yZXNvbHZlKHBhdGguZGlybmFtZShmaWxlUGF0aCksIHNyYylcbiAgICAgIH1cblxuICAgICAgaWYgKHJlbGF0aXZpemUgJiYgKGZpbGVQYXRoICE9PSB1bmRlZmluZWQgfHwgc2F2ZVBhdGggIT09IHVuZGVmaW5lZCkpIHtcbiAgICAgICAgY29uc3QgZnAgPSBzYXZlUGF0aCAhPT0gdW5kZWZpbmVkID8gc2F2ZVBhdGggOiBmaWxlUGF0aCFcbiAgICAgICAgc3JjID0gcGF0aC5yZWxhdGl2ZShwYXRoLmRpcm5hbWUoZnApLCBzcmMpXG4gICAgICB9XG5cbiAgICAgIC8vIFdhdGNoIGltYWdlIGZvciBjaGFuZ2VzXG4gICAgICBpZiAoaW1hZ2VXYXRjaGVyKSB7XG4gICAgICAgIGNvbnN0IHYgPSBpbWFnZVdhdGNoZXIud2F0Y2goc3JjKVxuICAgICAgICBpZiAodiAhPT0gdW5kZWZpbmVkKSBzcmMgPSBgJHtzcmN9P3Y9JHt2fWBcbiAgICAgIH1cblxuICAgICAgaW1nW2F0dHJOYW1lXSA9IHNyY1xuICAgIH1cbiAgfSlcbn1cblxuYXN5bmMgZnVuY3Rpb24gaGlnaGxpZ2h0Q29kZUJsb2NrcyhcbiAgZG9tRnJhZ21lbnQ6IERvY3VtZW50LFxuICBkZWZhdWx0TGFuZ3VhZ2U6IHN0cmluZyxcbikge1xuICBjb25zdCBmb250RmFtaWx5ID0gYXRvbS5jb25maWcuZ2V0KCdlZGl0b3IuZm9udEZhbWlseScpXG4gIGlmIChmb250RmFtaWx5KSB7XG4gICAgZm9yIChjb25zdCBjb2RlRWxlbWVudCBvZiBBcnJheS5mcm9tKFxuICAgICAgZG9tRnJhZ21lbnQucXVlcnlTZWxlY3RvckFsbCgnY29kZScpLFxuICAgICkpIHtcbiAgICAgIGNvZGVFbGVtZW50LnN0eWxlLmZvbnRGYW1pbHkgPSBmb250RmFtaWx5XG4gICAgfVxuICB9XG5cbiAgYXdhaXQgUHJvbWlzZS5hbGwoXG4gICAgQXJyYXkuZnJvbShkb21GcmFnbWVudC5xdWVyeVNlbGVjdG9yQWxsKCdwcmUnKSkubWFwKGFzeW5jIChwcmVFbGVtZW50KSA9PiB7XG4gICAgICBjb25zdCBjb2RlQmxvY2sgPVxuICAgICAgICBwcmVFbGVtZW50LmZpcnN0RWxlbWVudENoaWxkICE9PSBudWxsXG4gICAgICAgICAgPyBwcmVFbGVtZW50LmZpcnN0RWxlbWVudENoaWxkXG4gICAgICAgICAgOiBwcmVFbGVtZW50XG4gICAgICBjb25zdCBjYkNsYXNzID0gY29kZUJsb2NrLmNsYXNzTmFtZSB8fCBwcmVFbGVtZW50LmNsYXNzTmFtZVxuICAgICAgY29uc3QgZmVuY2VOYW1lID0gY2JDbGFzc1xuICAgICAgICA/IGNiQ2xhc3MucmVwbGFjZSgvXihsYW5nLXxzb3VyY2VDb2RlICkvLCAnJylcbiAgICAgICAgOiBkZWZhdWx0TGFuZ3VhZ2VcblxuICAgICAgY29uc3QgY3R3ID0gYXRvbUNvbmZpZygpLmNvZGVUYWJXaWR0aFxuICAgICAgY29uc3QgZWQgPSBuZXcgVGV4dEVkaXRvcih7XG4gICAgICAgIHJlYWRvbmx5OiB0cnVlLFxuICAgICAgICBrZXlib2FyZElucHV0RW5hYmxlZDogZmFsc2UsXG4gICAgICAgIHNob3dJbnZpc2libGVzOiBmYWxzZSxcbiAgICAgICAgdGFiTGVuZ3RoOiBjdHcgPT09IDAgPyBhdG9tLmNvbmZpZy5nZXQoJ2VkaXRvci50YWJMZW5ndGgnKSA6IGN0dyxcbiAgICAgIH0pXG4gICAgICBjb25zdCBlbCA9IGF0b20udmlld3MuZ2V0VmlldyhlZClcbiAgICAgIHRyeSB7XG4gICAgICAgIGVsLnNldFVwZGF0ZWRTeW5jaHJvbm91c2x5KHRydWUpXG4gICAgICAgIGVsLnN0eWxlLnBvaW50ZXJFdmVudHMgPSAnbm9uZSdcbiAgICAgICAgZWwuc3R5bGUucG9zaXRpb24gPSAnYWJzb2x1dGUnXG4gICAgICAgIGVsLnN0eWxlLnRvcCA9ICcxMDB2aCdcbiAgICAgICAgZWwuc3R5bGUud2lkdGggPSAnMTAwdncnXG4gICAgICAgIGF0b20uZ3JhbW1hcnMuYXNzaWduTGFuZ3VhZ2VNb2RlKFxuICAgICAgICAgIGVkLmdldEJ1ZmZlcigpLFxuICAgICAgICAgIHNjb3BlRm9yRmVuY2VOYW1lKGZlbmNlTmFtZSksXG4gICAgICAgIClcbiAgICAgICAgZWQuc2V0VGV4dChjb2RlQmxvY2sudGV4dENvbnRlbnQhLnJlcGxhY2UoL1xccj9cXG4kLywgJycpKVxuICAgICAgICBhdG9tLnZpZXdzLmdldFZpZXcoYXRvbS53b3Jrc3BhY2UpLmFwcGVuZENoaWxkKGVsKVxuICAgICAgICBhd2FpdCBlZGl0b3JUb2tlbml6ZWQoZWQpXG4gICAgICAgIGNvbnN0IGh0bWwgPSBBcnJheS5mcm9tKGVsLnF1ZXJ5U2VsZWN0b3JBbGwoJy5saW5lOm5vdCguZHVtbXkpJykpXG4gICAgICAgIHByZUVsZW1lbnQuY2xhc3NMaXN0LmFkZCgnZWRpdG9yLWNvbG9ycycpXG4gICAgICAgIHByZUVsZW1lbnQuaW5uZXJIVE1MID0gaHRtbC5tYXAoKHgpID0+IHguaW5uZXJIVE1MKS5qb2luKCdcXG4nKVxuICAgICAgICBpZiAoZmVuY2VOYW1lKSBwcmVFbGVtZW50LmNsYXNzTGlzdC5hZGQoYGxhbmctJHtmZW5jZU5hbWV9YClcbiAgICAgIH0gZmluYWxseSB7XG4gICAgICAgIGVsLnJlbW92ZSgpXG4gICAgICB9XG4gICAgfSksXG4gIClcblxuICByZXR1cm4gZG9tRnJhZ21lbnRcbn1cblxuYXN5bmMgZnVuY3Rpb24gZWRpdG9yVG9rZW5pemVkKGVkaXRvcjogVGV4dEVkaXRvcikge1xuICByZXR1cm4gbmV3IFByb21pc2UoKHJlc29sdmUpID0+IHtcbiAgICBjb25zdCBsYW5ndWFnZU1vZGUgPSBlZGl0b3IuZ2V0QnVmZmVyKCkuZ2V0TGFuZ3VhZ2VNb2RlKClcbiAgICBjb25zdCBuZXh0VXBkYXRlUHJvbWlzZSA9IGVkaXRvci5jb21wb25lbnQuZ2V0TmV4dFVwZGF0ZVByb21pc2UoKVxuICAgIGlmIChsYW5ndWFnZU1vZGUuZnVsbHlUb2tlbml6ZWQgfHwgbGFuZ3VhZ2VNb2RlLnRyZWUpIHtcbiAgICAgIHJlc29sdmUobmV4dFVwZGF0ZVByb21pc2UpXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbnN0IGRpc3AgPSBlZGl0b3Iub25EaWRUb2tlbml6ZSgoKSA9PiB7XG4gICAgICAgIGRpc3AuZGlzcG9zZSgpXG4gICAgICAgIHJlc29sdmUobmV4dFVwZGF0ZVByb21pc2UpXG4gICAgICB9KVxuICAgIH1cbiAgfSlcbn1cbiJdfQ==